[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$EngineBin,

    [Parameter(Mandatory = $true)]
    [string]$DwsimRevision,

    [Parameter(Mandatory = $true)]
    [string]$CasePath,

    [Parameter(Mandatory = $true)]
    [string]$OutputPath,

    [string]$PropertyPackageName,

    [string]$DirectPropertyPackageClass,

    [string]$DirectCompoundNamesCsv,

    [string]$CaseId,

    [Parameter(Mandatory = $true)]
    [double]$TemperatureK,

    [Parameter(Mandatory = $true)]
    [double]$PressurePa,

    [Parameter(Mandatory = $true)]
    [string]$CompositionCsv
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$Composition = @(
    $CompositionCsv.Split(",") | ForEach-Object {
        [double]::Parse($_, [Globalization.CultureInfo]::InvariantCulture)
    }
)

Add-Type -TypeDefinition @"
using System;
using System.Collections;
using System.Globalization;
using System.Linq;
using System.Reflection;

public static class DwsimThermoProbe
{
    private const BindingFlags PublicInstance =
        BindingFlags.Public | BindingFlags.Instance;

    public static object Get(object target, string name)
    {
        if (target == null) return null;
        PropertyInfo property = target.GetType().GetProperties(PublicInstance)
            .FirstOrDefault(candidate =>
                candidate.Name == name && candidate.GetIndexParameters().Length == 0);
        if (property != null) return property.GetValue(target, null);
        FieldInfo field = target.GetType().GetField(name, PublicInstance);
        return field == null ? null : field.GetValue(target);
    }

    public static object FirstValue(object dictionary)
    {
        if (!(dictionary is IEnumerable))
            throw new InvalidOperationException("Expected an enumerable dictionary.");
        foreach (object entry in (IEnumerable)dictionary)
        {
            object value = Get(entry, "Value");
            if (value != null) return value;
        }
        throw new InvalidOperationException("The dictionary is empty.");
    }

    private static object Invoke(object target, string name, object[] arguments, int count)
    {
        MethodInfo method = target.GetType().GetMethods(PublicInstance)
            .First(candidate =>
                candidate.Name == name && candidate.GetParameters().Length == count);
        return method.Invoke(target, arguments);
    }

    public static object CreateProbeFlowsheet(object automation, string[] compoundNames)
    {
        object flowsheet = Invoke(automation, "CreateFlowsheet", new object[0], 0);
        foreach (string compoundName in compoundNames)
            Invoke(flowsheet, "AddCompound", new object[] { compoundName }, 1);
        MethodInfo addObject = flowsheet.GetType().GetMethods(PublicInstance)
            .First(candidate => candidate.Name == "AddObject" && candidate.GetParameters().Length == 4);
        object materialStreamType = Enum.Parse(
            addObject.GetParameters()[0].ParameterType, "MaterialStream", true);
        addObject.Invoke(
            flowsheet,
            new object[] { materialStreamType, 50, 50, "thermo probe" });
        return flowsheet;
    }

    public static object CreatePropertyPackage(
        string fullTypeName, string componentName, object flowsheet)
    {
        Type packageType = AppDomain.CurrentDomain.GetAssemblies()
            .Select(assembly => assembly.GetType(fullTypeName, false, false))
            .FirstOrDefault(candidate => candidate != null);
        if (packageType == null)
            throw new InvalidOperationException(
                "The property-package class is not loaded: " + fullTypeName);
        object propertyPackage = Activator.CreateInstance(packageType);
        Set(propertyPackage, "ComponentName", componentName);
        Set(propertyPackage, "Flowsheet", flowsheet);
        return propertyPackage;
    }

    private static void Set(object target, string name, object value)
    {
        PropertyInfo property = target.GetType().GetProperties(PublicInstance)
            .First(candidate =>
                candidate.Name == name && candidate.CanWrite &&
                candidate.GetIndexParameters().Length == 0);
        property.SetValue(target, value, null);
    }

    public static object ValueByComponentName(object dictionary, string componentName)
    {
        if (!(dictionary is IEnumerable))
            throw new InvalidOperationException("Expected an enumerable dictionary.");
        foreach (object entry in (IEnumerable)dictionary)
        {
            object value = Get(entry, "Value");
            string candidate = Convert.ToString(
                Get(value, "ComponentName"), CultureInfo.InvariantCulture);
            if (value != null && String.Equals(
                    candidate, componentName, StringComparison.Ordinal)) return value;
        }
        throw new InvalidOperationException(
            "The property-package dictionary has no exact ComponentName match for '" +
            componentName + "'.");
    }

    public static string FullTypeName(object target)
    {
        return target == null ? null : target.GetType().FullName;
    }

    public static object FirstMaterialStream(object flowsheet)
    {
        object objects = Get(flowsheet, "SimulationObjects");
        if (!(objects is IEnumerable))
            throw new InvalidOperationException("Flowsheet objects are not enumerable.");
        foreach (object entry in (IEnumerable)objects)
        {
            object value = Get(entry, "Value");
            if (value != null && value.GetType().FullName ==
                "DWSIM.Thermodynamics.Streams.MaterialStream") return value;
        }
        throw new InvalidOperationException("The flowsheet has no material stream.");
    }

    public static void SetCurrentStream(object propertyPackage, object stream)
    {
        PropertyInfo property = propertyPackage.GetType().GetProperties(PublicInstance)
            .First(candidate => candidate.Name == "CurrentMaterialStream" && candidate.CanWrite);
        property.SetValue(propertyPackage, stream, null);
    }

    public static string[] CompoundNames(object propertyPackage)
    {
        MethodInfo method = propertyPackage.GetType().GetMethods(PublicInstance)
            .First(candidate => candidate.Name == "RET_VNAMES" && candidate.GetParameters().Length == 0);
        return ((IEnumerable)method.Invoke(propertyPackage, null))
            .Cast<object>()
            .Select(value => Convert.ToString(value, CultureInfo.InvariantCulture))
            .ToArray();
    }

    public static double[] FugacityCoefficients(
        object propertyPackage, double[] composition, double temperature, double pressure,
        string stateName)
    {
        MethodInfo method = propertyPackage.GetType().GetMethods(PublicInstance)
            .First(candidate =>
                candidate.Name == "DW_CalcFugCoeff" &&
                candidate.GetParameters().Length == 4 &&
                candidate.GetParameters()[0].ParameterType == typeof(Array));
        object state = Enum.Parse(method.GetParameters()[3].ParameterType, stateName, true);
        return ((IEnumerable)method.Invoke(
                propertyPackage,
                new object[] { composition, temperature, pressure, state }))
            .Cast<object>()
            .Select(value => Convert.ToDouble(value, CultureInfo.InvariantCulture))
            .ToArray();
    }

    public static object[] FlashPT(
        object propertyPackage, double[] composition, double temperature, double pressure)
    {
        object flash = Get(propertyPackage, "FlashBase");
        MethodInfo method = flash.GetType().GetMethods(PublicInstance)
            .First(candidate =>
                candidate.Name == "Flash_PT" && candidate.GetParameters().Length == 6);
        return ((IEnumerable)method.Invoke(
                flash,
                new object[] { composition, pressure, temperature, propertyPackage, false, null }))
            .Cast<object>()
            .ToArray();
    }

    public static double[] DoubleArray(object value)
    {
        return ((IEnumerable)value).Cast<object>()
            .Select(item => Convert.ToDouble(item, CultureInfo.InvariantCulture))
            .ToArray();
    }
}
"@

$engineDirectory = (Resolve-Path -LiteralPath $EngineBin).Path
$resolvedCase = (Resolve-Path -LiteralPath $CasePath).Path
$resolvedCaseId = if ([string]::IsNullOrWhiteSpace($CaseId)) {
    [IO.Path]::GetFileNameWithoutExtension($resolvedCase)
}
else {
    $CaseId
}
$automationPath = Join-Path $engineDirectory "DWSIM.Automation.dll"
$interfacesPath = Join-Path $engineDirectory "DWSIM.Interfaces.dll"

if (-not (Test-Path -LiteralPath $automationPath)) {
    throw "Missing DWSIM.Automation.dll in $engineDirectory"
}
if (-not (Test-Path -LiteralPath $interfacesPath)) {
    throw "Missing DWSIM.Interfaces.dll in $engineDirectory"
}
if ([string]::IsNullOrWhiteSpace($DwsimRevision)) {
    throw "DwsimRevision is required for an auditable capture"
}

[Reflection.Assembly]::LoadFrom($interfacesPath) | Out-Null
$thermodynamicsPath = Join-Path $engineDirectory "DWSIM.Thermodynamics.dll"
if (-not (Test-Path -LiteralPath $thermodynamicsPath)) {
    throw "Missing DWSIM.Thermodynamics.dll in $engineDirectory"
}
[Reflection.Assembly]::LoadFrom($thermodynamicsPath) | Out-Null
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

$automation = New-Object DWSIM.Automation.Automation3
$directConstruction = -not [string]::IsNullOrWhiteSpace($DirectPropertyPackageClass)
if ($directConstruction) {
    if ([string]::IsNullOrWhiteSpace($PropertyPackageName)) {
        throw "PropertyPackageName is required with DirectPropertyPackageClass"
    }
    if ([string]::IsNullOrWhiteSpace($DirectCompoundNamesCsv)) {
        throw "DirectCompoundNamesCsv is required with DirectPropertyPackageClass"
    }
    $directCompoundNames = @(
        $DirectCompoundNamesCsv.Split(",") | ForEach-Object { $_.Trim() }
    )
    if ($directCompoundNames.Count -eq 0 -or $directCompoundNames -contains "") {
        throw "DirectCompoundNamesCsv must contain non-empty DWSIM catalog names"
    }
    $flowsheet = [DwsimThermoProbe]::CreateProbeFlowsheet(
        $automation, $directCompoundNames
    )
    $propertyPackage = [DwsimThermoProbe]::CreatePropertyPackage(
        $DirectPropertyPackageClass, $PropertyPackageName, $flowsheet
    )
}
else {
    $flowsheet = $automation.LoadFlowsheet2($resolvedCase)
    $packages = [DwsimThermoProbe]::Get($flowsheet, "PropertyPackages")
    $propertyPackage = if ([string]::IsNullOrWhiteSpace($PropertyPackageName)) {
        [DwsimThermoProbe]::FirstValue($packages)
    }
    else {
        [DwsimThermoProbe]::ValueByComponentName($packages, $PropertyPackageName)
    }
}
$materialStream = [DwsimThermoProbe]::FirstMaterialStream($flowsheet)
[DwsimThermoProbe]::SetCurrentStream($propertyPackage, $materialStream)

$compoundNames = @([DwsimThermoProbe]::CompoundNames($propertyPackage))
if ($Composition.Count -ne $compoundNames.Count) {
    throw "Composition length must match the DWSIM package compound count"
}
if ([Math]::Abs(($Composition | Measure-Object -Sum).Sum - 1.0) -gt 1.0e-12) {
    throw "Composition must sum to one"
}

$liquidFugacity = @([DwsimThermoProbe]::FugacityCoefficients(
    $propertyPackage, $Composition, $TemperatureK, $PressurePa, "Liquid"
))
$vaporFugacity = @([DwsimThermoProbe]::FugacityCoefficients(
    $propertyPackage, $Composition, $TemperatureK, $PressurePa, "Vapor"
))
$flash = [DwsimThermoProbe]::FlashPT(
    $propertyPackage, $Composition, $TemperatureK, $PressurePa
)

$document = [ordered]@{
    schema_version = "dwsim-thermo-package-golden-1"
    case_id = $resolvedCaseId
    source = [ordered]@{
        dwsim_revision = $DwsimRevision
        automation_version = [string]$automation.GetVersion()
        input_file_sha256 = (Get-FileHash -Algorithm SHA256 -Path $resolvedCase).Hash.ToLowerInvariant()
        property_package = [string][DwsimThermoProbe]::Get($propertyPackage, "ComponentName")
        property_package_class = [DwsimThermoProbe]::FullTypeName($propertyPackage)
        property_package_construction = if ($directConstruction) {
            "direct-class-over-case-compound-domain"
        }
        else {
            "deserialized-from-case"
        }
        flash_algorithm = [string][DwsimThermoProbe]::Get(
            [DwsimThermoProbe]::Get($propertyPackage, "FlashBase"), "Name"
        )
    }
    inputs = [ordered]@{
        compounds = $compoundNames
        temperature_k = $TemperatureK
        pressure_pa = $PressurePa
        composition = $Composition
    }
    outputs = [ordered]@{
        liquid_fugacity_coefficients = $liquidFugacity
        vapor_fugacity_coefficients = $vaporFugacity
        liquid_fraction = [double]$flash[0]
        vapor_fraction = [double]$flash[1]
        liquid_composition = @([DwsimThermoProbe]::DoubleArray($flash[2]))
        vapor_composition = @([DwsimThermoProbe]::DoubleArray($flash[3]))
        iterations = [int]$flash[4]
        equilibrium_ratios = @([DwsimThermoProbe]::DoubleArray($flash[9]))
    }
}

$outputDirectory = Split-Path -Parent $OutputPath
if (-not [string]::IsNullOrWhiteSpace($outputDirectory)) {
    New-Item -ItemType Directory -Path $outputDirectory -Force | Out-Null
}
$document | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $OutputPath -Encoding utf8
