[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$EngineBin,

    [Parameter(Mandatory = $true)]
    [string]$SourceCasePath,

    [Parameter(Mandatory = $true)]
    [string]$OutputCasePath,

    [Parameter(Mandatory = $true)]
    [string]$ObjectTag,

    [Parameter(Mandatory = $true)]
    [string]$PropertyName,

    [Parameter(Mandatory = $true)]
    [double]$PropertyValue
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

Add-Type -TypeDefinition @"
using System;
using System.Collections;
using System.Reflection;

public static class DwsimCaseMutationReflection
{
    private const BindingFlags PublicInstance =
        BindingFlags.Public | BindingFlags.Instance;

    public static object Get(object target, string name)
    {
        if (target == null) return null;
        PropertyInfo property = target.GetType().GetProperty(name, PublicInstance);
        return property == null ? null : property.GetValue(target, null);
    }

    public static object Invoke(object target, string name, object[] arguments)
    {
        if (target == null) throw new ArgumentNullException("target");
        return target.GetType().InvokeMember(
            name,
            PublicInstance | BindingFlags.InvokeMethod,
            null,
            target,
            arguments);
    }

    public static object FindByTag(object flowsheet, string requestedTag)
    {
        object simulationObjects = Get(flowsheet, "SimulationObjects");
        object values = Get(simulationObjects, "Values");
        if (!(values is IEnumerable))
            throw new InvalidOperationException("DWSIM flowsheet has no enumerable simulation objects.");

        foreach (object simulationObject in (IEnumerable)values)
        {
            string name = Convert.ToString(Get(simulationObject, "Name"));
            object graphicObject = Get(simulationObject, "GraphicObject");
            string tag = Convert.ToString(Get(graphicObject, "Tag"));
            if (String.Equals(requestedTag, tag, StringComparison.Ordinal) ||
                String.Equals(requestedTag, name, StringComparison.Ordinal))
                return simulationObject;
        }

        throw new InvalidOperationException("DWSIM object tag was not found: " + requestedTag);
    }
}
"@

$engineDirectory = (Resolve-Path $EngineBin).Path
$sourcePath = (Resolve-Path $SourceCasePath).Path
$outputPath = [IO.Path]::GetFullPath((Join-Path (Get-Location) $OutputCasePath))

if ([string]::Equals($sourcePath, $outputPath, [StringComparison]::OrdinalIgnoreCase)) {
    throw "Source and output DWSIM case paths must differ"
}

if (Test-Path -LiteralPath $outputPath) {
    throw "Output DWSIM case already exists: $outputPath"
}

$outputDirectory = Split-Path -Parent $outputPath
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
$thermoCAssemblyPath = Get-ChildItem -Path $engineDirectory -Filter "ThermoCS.dll" -File -Recurse -ErrorAction SilentlyContinue |
    Select-Object -First 1
if ($null -ne $thermoCAssemblyPath) {
    [Reflection.Assembly]::LoadFrom($thermoCAssemblyPath.FullName) | Out-Null
}
[Reflection.Assembly]::LoadFrom($automationPath) | Out-Null

$automation = New-Object DWSIM.Automation.Automation3
$flowsheet = $automation.LoadFlowsheet2($sourcePath)
$simulationObject = [DwsimCaseMutationReflection]::FindByTag($flowsheet, $ObjectTag)

$calculationType = [DwsimCaseMutationReflection]::Invoke(
    $simulationObject,
    "GetPropertyValue",
    [object[]]@("ThermalProfile,CalculationType", $null)
)
$heatBalanceEnabled = [bool][DwsimCaseMutationReflection]::Get(
    $simulationObject,
    "CalculateHeatBalance"
)

if (-not $heatBalanceEnabled) {
    throw "Target pipe does not have its heat balance enabled"
}
if ([string]$calculationType -ne "Definir_CGTC") {
    throw "Target pipe is not in defined-HTC mode: $calculationType"
}

$accepted = [bool][DwsimCaseMutationReflection]::Invoke(
    $simulationObject,
    "SetPropertyValue",
    [object[]]@($PropertyName, $PropertyValue, $null)
)
if (-not $accepted) {
    throw "DWSIM rejected property '$PropertyName'"
}

$errors = @($automation.CalculateFlowsheet4($flowsheet))
if ($errors.Count -ne 0) {
    $messages = @($errors | ForEach-Object { [string]$_.Message })
    throw "DWSIM solve failed: $($messages -join '; ')"
}

$automation.SaveFlowsheet2($flowsheet, $outputPath)

$savedFlowsheet = $automation.LoadFlowsheet2($outputPath)
$savedObject = [DwsimCaseMutationReflection]::FindByTag($savedFlowsheet, $ObjectTag)
$savedValue = [double][DwsimCaseMutationReflection]::Invoke(
    $savedObject,
    "GetPropertyValue",
    [object[]]@($PropertyName, $null)
)
$savedCalculationType = [DwsimCaseMutationReflection]::Invoke(
    $savedObject,
    "GetPropertyValue",
    [object[]]@("ThermalProfile,CalculationType", $null)
)
$savedHeatBalanceEnabled = [bool][DwsimCaseMutationReflection]::Get(
    $savedObject,
    "CalculateHeatBalance"
)

if ([Math]::Abs($savedValue - $PropertyValue) -gt 1e-12) {
    throw "Saved DWSIM property differs: expected $PropertyValue, found $savedValue"
}
if (-not $savedHeatBalanceEnabled -or [string]$savedCalculationType -ne "Definir_CGTC") {
    throw "Saved DWSIM thermal mode fields differ from the source case"
}

[ordered]@{
    source_case = $sourcePath
    output_case = $outputPath
    object_tag = $ObjectTag
    property = $PropertyName
    value = $savedValue
    heat_balance_enabled = $savedHeatBalanceEnabled
    calculation_type = [string]$savedCalculationType
    solve_error_count = $errors.Count
} | ConvertTo-Json -Depth 4
