[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$EngineBin,

    [Parameter(Mandatory = $true)]
    [string]$InputPath,

    [Parameter(Mandatory = $true)]
    [string]$OutputPath,

    [Parameter(Mandatory = $true)]
    [string]$PropertyPackageName,

    [Parameter(Mandatory = $true)]
    [string]$Pair,

    [Parameter(Mandatory = $true)]
    [AllowEmptyString()]
    [string]$Expression
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

Add-Type -TypeDefinition @"
using System;
using System.Collections;
using System.Linq;
using System.Reflection;

public static class DwsimAdvancedInteractionConfigurator
{
    private const BindingFlags PublicInstance = BindingFlags.Public | BindingFlags.Instance;

    public static object Get(object target, string name)
    {
        PropertyInfo property = target.GetType().GetProperties(PublicInstance)
            .FirstOrDefault(candidate =>
                candidate.Name == name && candidate.GetIndexParameters().Length == 0);
        if (property != null) return property.GetValue(target, null);
        FieldInfo field = target.GetType().GetField(name, PublicInstance);
        if (field != null) return field.GetValue(target);
        throw new InvalidOperationException("Missing reflected member: " + name);
    }

    public static object PackageByComponentName(object dictionary, string componentName)
    {
        foreach (object entry in (IEnumerable)dictionary)
        {
            object value = Get(entry, "Value");
            if (String.Equals(
                    Convert.ToString(Get(value, "ComponentName")),
                    componentName,
                    StringComparison.Ordinal)) return value;
        }
        throw new InvalidOperationException(
            "No exact property-package ComponentName match for '" + componentName + "'.");
    }

    public static void SetExpression(object propertyPackage, string pair, string expression)
    {
        object value = Get(propertyPackage, "KijExpressions");
        IDictionary dictionary = value as IDictionary;
        if (dictionary == null)
            throw new InvalidOperationException("KijExpressions is not an IDictionary.");
        dictionary[pair] = expression;
    }
}
"@

$engineDirectory = (Resolve-Path -LiteralPath $EngineBin).Path
$resolvedInput = (Resolve-Path -LiteralPath $InputPath).Path
$resolvedOutput = [IO.Path]::GetFullPath((Join-Path (Get-Location) $OutputPath))
if (Test-Path -LiteralPath $resolvedOutput) {
    throw "Output DWSIM case already exists: $resolvedOutput"
}
$outputDirectory = Split-Path -Parent $resolvedOutput
if (-not (Test-Path -LiteralPath $outputDirectory -PathType Container)) {
    throw "Output directory does not exist: $outputDirectory"
}
if ($Pair -notmatch "^[^/]+/[^/]+$") {
    throw "Pair must use exact DWSIM names in first/second form"
}

$interfacesPath = Join-Path $engineDirectory "DWSIM.Interfaces.dll"
$automationPath = Join-Path $engineDirectory "DWSIM.Automation.dll"
foreach ($assemblyPath in @($interfacesPath, $automationPath)) {
    if (-not (Test-Path -LiteralPath $assemblyPath -PathType Leaf)) {
        throw "Missing DWSIM assembly: $assemblyPath"
    }
    [Reflection.Assembly]::LoadFrom($assemblyPath) | Out-Null
}
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

$automation = New-Object DWSIM.Automation.Automation3
$flowsheet = $automation.LoadFlowsheet2($resolvedInput)
$packages = [DwsimAdvancedInteractionConfigurator]::Get(
    $flowsheet, "PropertyPackages"
)
$propertyPackage = [DwsimAdvancedInteractionConfigurator]::PackageByComponentName(
    $packages, $PropertyPackageName
)
[DwsimAdvancedInteractionConfigurator]::SetExpression(
    $propertyPackage, $Pair, $Expression
)
$automation.SaveFlowsheet2($flowsheet, $resolvedOutput)

[ordered]@{
    output_case = $resolvedOutput
    property_package = $PropertyPackageName
    pair = $Pair
    expression = $Expression
    sha256 = (Get-FileHash -LiteralPath $resolvedOutput -Algorithm SHA256).Hash.ToLowerInvariant()
} | ConvertTo-Json
