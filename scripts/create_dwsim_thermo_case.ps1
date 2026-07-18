[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$EngineBin,

    [Parameter(Mandatory = $true)]
    [string]$OutputPath,

    [Parameter(Mandatory = $true)]
    [string]$PropertyPackageName,

    [Parameter(Mandatory = $true)]
    [string]$CompoundNamesCsv,

    [string]$DirectPropertyPackageClass,

    [string]$CompoundPropertyOverridesPath
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

Add-Type -TypeDefinition @"
using System;
using System.Collections;
using System.Globalization;
using System.Linq;
using System.Reflection;

public static class DwsimThermoCaseBuilder
{
    private const BindingFlags PublicInstance = BindingFlags.Public | BindingFlags.Instance;

    private static object Invoke(object target, string name, object[] arguments, int count)
    {
        MethodInfo method = target.GetType().GetMethods(PublicInstance)
            .First(candidate => candidate.Name == name && candidate.GetParameters().Length == count);
        return method.Invoke(target, arguments);
    }

    public static object Create(
        object automation, string[] compoundNames, string propertyPackageName,
        string directPropertyPackageClass)
    {
        object flowsheet = Invoke(automation, "CreateFlowsheet", new object[0], 0);
        foreach (string compoundName in compoundNames)
            Invoke(flowsheet, "AddCompound", new object[] { compoundName }, 1);

        MethodInfo addObject = flowsheet.GetType().GetMethods(PublicInstance)
            .First(candidate => candidate.Name == "AddObject" && candidate.GetParameters().Length == 4);
        object materialStreamType = Enum.Parse(
            addObject.GetParameters()[0].ParameterType, "MaterialStream", true);
        Invoke(
            flowsheet,
            "AddObject",
            new object[] { materialStreamType, 50, 50, "thermo probe" },
            4);
        if (String.IsNullOrWhiteSpace(directPropertyPackageClass))
            Invoke(flowsheet,"CreateAndAddPropertyPackage",new object[] { propertyPackageName },1);
        else
        {
            Type packageType=AppDomain.CurrentDomain.GetAssemblies()
                .Select(assembly=>assembly.GetType(directPropertyPackageClass,false,false))
                .FirstOrDefault(candidate=>candidate!=null);
            if(packageType==null)throw new InvalidOperationException(
                "Property-package class is not loaded: "+directPropertyPackageClass);
            object package=Activator.CreateInstance(packageType);
            packageType.GetProperties(PublicInstance)
                .First(property=>property.Name=="ComponentName"&&property.CanWrite&&property.GetIndexParameters().Length==0)
                .SetValue(package,propertyPackageName,null);
            Invoke(flowsheet,"AddPropertyPackage",new object[] { package },1);
        }
        return flowsheet;
    }

    public static void SetCompoundProperty(
        object flowsheet, string compoundName, string propertyName, double value)
    {
        object selected = flowsheet.GetType().GetProperties(PublicInstance)
            .First(property => property.Name == "SelectedCompounds" &&
                property.GetIndexParameters().Length == 0)
            .GetValue(flowsheet, null);
        foreach (object entry in (IEnumerable)selected)
        {
            object constantProperties = entry.GetType().GetProperty("Value")
                .GetValue(entry, null);
            string name = Convert.ToString(
                constantProperties.GetType().GetProperty("Name")
                    .GetValue(constantProperties, null),
                CultureInfo.InvariantCulture);
            if (!String.Equals(name, compoundName, StringComparison.Ordinal)) continue;
            PropertyInfo property = constantProperties.GetType().GetProperties(PublicInstance)
                .First(candidate => candidate.Name == propertyName && candidate.CanWrite &&
                    candidate.GetIndexParameters().Length == 0);
            property.SetValue(constantProperties, value, null);
            return;
        }
        throw new InvalidOperationException(
            "Selected compound was not found for override: " + compoundName);
    }
}
"@

$engineDirectory = (Resolve-Path -LiteralPath $EngineBin).Path
$outputFile = [IO.Path]::GetFullPath((Join-Path (Get-Location) $OutputPath))
if (Test-Path -LiteralPath $outputFile) {
    throw "Output DWSIM case already exists: $outputFile"
}
$outputDirectory = Split-Path -Parent $outputFile
if (-not (Test-Path -LiteralPath $outputDirectory -PathType Container)) {
    throw "Output directory does not exist: $outputDirectory"
}

$interfacesPath = Join-Path $engineDirectory "DWSIM.Interfaces.dll"
$automationPath = Join-Path $engineDirectory "DWSIM.Automation.dll"
$thermodynamicsPath = Join-Path $engineDirectory "DWSIM.Thermodynamics.dll"
foreach ($assemblyPath in @($interfacesPath, $automationPath, $thermodynamicsPath)) {
    if (-not (Test-Path -LiteralPath $assemblyPath -PathType Leaf)) {
        throw "Missing DWSIM assembly: $assemblyPath"
    }
}
[Reflection.Assembly]::LoadFrom($interfacesPath) | Out-Null
$thermoCAssemblyPath = Get-ChildItem `
    -Path $engineDirectory `
    -Filter "ThermoCS.dll" `
    -File `
    -Recurse `
    -ErrorAction SilentlyContinue |
    Select-Object -First 1
if ($null -ne $thermoCAssemblyPath) {
    [Reflection.Assembly]::LoadFrom($thermoCAssemblyPath.FullName) | Out-Null
}
[Reflection.Assembly]::LoadFrom($thermodynamicsPath) | Out-Null
[Reflection.Assembly]::LoadFrom($automationPath) | Out-Null

$compoundNames = @(
    $CompoundNamesCsv.Split(",") | ForEach-Object { $_.Trim() }
)
if ($compoundNames.Count -eq 0 -or $compoundNames -contains "") {
    throw "CompoundNamesCsv must contain non-empty DWSIM catalog names"
}

$automation = New-Object DWSIM.Automation.Automation3
$flowsheet = [DwsimThermoCaseBuilder]::Create(
    $automation, $compoundNames, $PropertyPackageName, $DirectPropertyPackageClass
)
if (-not [string]::IsNullOrWhiteSpace($CompoundPropertyOverridesPath)) {
    $overrideFile = (Resolve-Path -LiteralPath $CompoundPropertyOverridesPath).Path
    $overrides = Get-Content -LiteralPath $overrideFile -Raw | ConvertFrom-Json
    foreach ($compoundProperty in $overrides.PSObject.Properties) {
        if ($compoundNames -cnotcontains $compoundProperty.Name) {
            throw "Override compound is outside the exact case domain: $($compoundProperty.Name)"
        }
        foreach ($property in $compoundProperty.Value.PSObject.Properties) {
            [DwsimThermoCaseBuilder]::SetCompoundProperty(
                $flowsheet, $compoundProperty.Name, $property.Name,
                [double]$property.Value
            )
        }
    }
}
$automation.SaveFlowsheet2($flowsheet, $outputFile)

[ordered]@{
    output_case = $outputFile
    property_package = $PropertyPackageName
    compounds = $compoundNames
    compound_property_overrides = if ([string]::IsNullOrWhiteSpace($CompoundPropertyOverridesPath)) { $null } else { (Resolve-Path -LiteralPath $CompoundPropertyOverridesPath).Path }
    sha256 = (Get-FileHash -Algorithm SHA256 -Path $outputFile).Hash.ToLowerInvariant()
} | ConvertTo-Json -Depth 4
