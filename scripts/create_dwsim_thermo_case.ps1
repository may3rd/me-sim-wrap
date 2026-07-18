[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$EngineBin,

    [Parameter(Mandatory = $true)]
    [string]$OutputPath,

    [Parameter(Mandatory = $true)]
    [string]$PropertyPackageName,

    [Parameter(Mandatory = $true)]
    [string]$CompoundNamesCsv
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

Add-Type -TypeDefinition @"
using System;
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
        object automation, string[] compoundNames, string propertyPackageName)
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
        Invoke(
            flowsheet,
            "CreateAndAddPropertyPackage",
            new object[] { propertyPackageName },
            1);
        return flowsheet;
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
foreach ($assemblyPath in @($interfacesPath, $automationPath)) {
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
[Reflection.Assembly]::LoadFrom($automationPath) | Out-Null

$compoundNames = @(
    $CompoundNamesCsv.Split(",") | ForEach-Object { $_.Trim() }
)
if ($compoundNames.Count -eq 0 -or $compoundNames -contains "") {
    throw "CompoundNamesCsv must contain non-empty DWSIM catalog names"
}

$automation = New-Object DWSIM.Automation.Automation3
$flowsheet = [DwsimThermoCaseBuilder]::Create(
    $automation, $compoundNames, $PropertyPackageName
)
$automation.SaveFlowsheet2($flowsheet, $outputFile)

[ordered]@{
    output_case = $outputFile
    property_package = $PropertyPackageName
    compounds = $compoundNames
    sha256 = (Get-FileHash -Algorithm SHA256 -Path $outputFile).Hash.ToLowerInvariant()
} | ConvertTo-Json -Depth 4
