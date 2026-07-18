[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$EngineBin,

    [Parameter(Mandatory = $true)]
    [string]$DwsimRevision,

    [Parameter(Mandatory = $true)]
    [string]$OutputPath,

    [string]$CasePath,
    [string]$CaseId,
    [string]$PropertyPackage = "unknown",
    [string]$FlashAlgorithm = "unknown",
    [switch]$CalculateBubbleAndDewPoints,
    [switch]$CaptureColumnProfiles,
    [string[]]$ObjectTags = @(),
    [string[]]$PropertyNames = @(),
    [string[]]$CompoundNames = @(
        "Methane",
        "Ethane",
        "Propane",
        "N-butane",
        "N-pentane"
    )
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

Add-Type -TypeDefinition @"
using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Reflection;

public static class DwsimCaptureReflection
{
    private const BindingFlags PublicInstance =
        BindingFlags.Public | BindingFlags.Instance;

    public static object Get(object target, string name)
    {
        if (target == null)
        {
            return null;
        }

        PropertyInfo property = target.GetType().GetProperties(PublicInstance)
            .FirstOrDefault(candidate =>
                candidate.Name == name &&
                candidate.GetIndexParameters().Length == 0);

        if (property != null)
        {
            return property.GetValue(target, null);
        }

        FieldInfo field = target.GetType().GetField(name, PublicInstance);
        return field == null ? null : field.GetValue(target);
    }

    public static object Invoke(
        object target,
        string name,
        object[] arguments)
    {
        if (target == null)
        {
            throw new ArgumentNullException("target");
        }

        return target.GetType().InvokeMember(
            name,
            PublicInstance | BindingFlags.InvokeMethod,
            null,
            target,
            arguments);
    }

    public static string TypeName(object target)
    {
        return target == null
            ? null
            : target.GetType().Name;
    }

    public static IDictionary StringDoubleDictionary(object target)
    {
        SortedDictionary<string, double> result =
            new SortedDictionary<string, double>(StringComparer.Ordinal);

        if (!(target is IEnumerable))
        {
            return result;
        }

        foreach (object entry in (IEnumerable)target)
        {
            object key = Get(entry, "Key");
            object value = Get(entry, "Value");

            if (key != null && value != null)
            {
                result[Convert.ToString(key, CultureInfo.InvariantCulture)] =
                    Convert.ToDouble(value, CultureInfo.InvariantCulture);
            }
        }

        return result;
    }

    public static object NumericArray(object target)
    {
        if (target == null || target is string)
        {
            return null;
        }

        if (target is IEnumerable)
        {
            List<object> result = new List<object>();
            foreach (object item in (IEnumerable)target)
            {
                result.Add(item is IEnumerable && !(item is string)
                    ? NumericArray(item)
                    : (object)Convert.ToDouble(item, CultureInfo.InvariantCulture));
            }
            return result.ToArray();
        }

        return null;
    }

    public static object ColumnEnthalpyArray(
        object target,
        string compositionField,
        string stateName)
    {
        object propertyPackage = Get(target, "PropertyPackage");
        object temperatures = Get(target, "Tf");
        object pressures = Get(target, "P0");
        object compositions = Get(target, compositionField);
        if (
            propertyPackage == null ||
            !(temperatures is IEnumerable) ||
            !(pressures is IEnumerable) ||
            !(compositions is IEnumerable))
        {
            return null;
        }

        MethodInfo method = propertyPackage.GetType().GetMethods(PublicInstance)
            .FirstOrDefault(candidate =>
                candidate.Name == "DW_CalcEnthalpy" &&
                candidate.GetParameters().Length == 4 &&
                candidate.GetParameters()[0].ParameterType == typeof(Array));
        if (method == null)
        {
            return null;
        }
        Type stateType = method.GetParameters()[3].ParameterType;
        object state = Enum.Parse(stateType, stateName, true);
        PropertyInfo currentStreamProperty = propertyPackage.GetType()
            .GetProperties(PublicInstance)
            .FirstOrDefault(candidate =>
                candidate.Name == "CurrentMaterialStream" &&
                candidate.CanRead &&
                candidate.CanWrite);
        object originalStream = currentStreamProperty == null
            ? null
            : currentStreamProperty.GetValue(propertyPackage, null);
        try
        {
            if (currentStreamProperty != null && originalStream == null)
            {
                object flowsheet = Get(target, "FlowSheet");
                object simulationObjects = Get(flowsheet, "SimulationObjects");
                if (simulationObjects is IEnumerable)
                {
                    foreach (object entry in (IEnumerable)simulationObjects)
                    {
                        object candidate = Get(entry, "Value");
                        if (
                            candidate != null &&
                            candidate.GetType().FullName ==
                                "DWSIM.Thermodynamics.Streams.MaterialStream")
                        {
                            currentStreamProperty.SetValue(
                                propertyPackage,
                                candidate,
                                null);
                            break;
                        }
                    }
                }
            }
            List<double> temperatureValues = ((IEnumerable)temperatures)
                .Cast<object>()
                .Select(item => Convert.ToDouble(item, CultureInfo.InvariantCulture))
                .ToList();
            List<double> pressureValues = ((IEnumerable)pressures)
                .Cast<object>()
                .Select(item => Convert.ToDouble(item, CultureInfo.InvariantCulture))
                .ToList();
            List<object> compositionValues = ((IEnumerable)compositions)
                .Cast<object>()
                .ToList();
            if (
                temperatureValues.Count != pressureValues.Count ||
                temperatureValues.Count != compositionValues.Count)
            {
                return null;
            }

            object[] result = new object[temperatureValues.Count];
            for (int index = 0; index < result.Length; index++)
            {
                double[] composition = ((IEnumerable)compositionValues[index])
                    .Cast<object>()
                    .Select(item => Convert.ToDouble(item, CultureInfo.InvariantCulture))
                    .ToArray();
                result[index] = Convert.ToDouble(
                    method.Invoke(
                        propertyPackage,
                        new object[] {
                            composition,
                            temperatureValues[index],
                            pressureValues[index],
                            state
                        }),
                    CultureInfo.InvariantCulture);
            }
            return result;
        }
        finally
        {
            if (currentStreamProperty != null)
            {
                currentStreamProperty.SetValue(propertyPackage, originalStream, null);
            }
        }
    }

    public static int SetFlashSetting(
        object flowsheet,
        string settingName,
        string settingValue)
    {
        object packages = Get(flowsheet, "PropertyPackages");
        object values = Get(packages, "Values");
        int count = 0;

        if (!(values is IEnumerable))
        {
            throw new InvalidOperationException(
                "DWSIM flowsheet has no enumerable property packages.");
        }

        foreach (object package in (IEnumerable)values)
        {
            object settings = Get(package, "FlashSettings");
            Type[] arguments = settings == null
                ? new Type[0]
                : settings.GetType().GetGenericArguments();
            PropertyInfo item = settings == null
                ? null
                : settings.GetType().GetProperty("Item", PublicInstance);

            if (arguments.Length != 2 || !arguments[0].IsEnum || item == null)
            {
                throw new InvalidOperationException(
                    "DWSIM property package has no writable FlashSettings dictionary.");
            }

            object key = Enum.Parse(arguments[0], settingName);
            item.SetValue(settings, settingValue, new object[] { key });
            count++;
        }

        if (count == 0)
        {
            throw new InvalidOperationException(
                "DWSIM flowsheet has no property packages.");
        }

        return count;
    }
}
"@

function Get-MemberValue {
    param(
        [AllowNull()]
        [object]$Object,

        [Parameter(Mandatory = $true)]
        [string]$Name
    )

    if ($null -eq $Object) {
        return $null
    }

    $property = $Object.PSObject.Properties[$Name]

    if ($null -eq $property) {
        return $null
    }

    return $property.Value
}

function Get-NumericText {
    param(
        [Parameter(Mandatory = $true)]
        [object]$Value
    )

    $culture = [Globalization.CultureInfo]::InvariantCulture

    if ($Value -is [double]) {
        return $Value.ToString("R", $culture)
    }

    if ($Value -is [float]) {
        return $Value.ToString("R", $culture)
    }

    if ($Value -is [decimal]) {
        return $Value.ToString("G29", $culture)
    }

    if (
        $Value -is [byte] -or
        $Value -is [sbyte] -or
        $Value -is [int16] -or
        $Value -is [uint16] -or
        $Value -is [int32] -or
        $Value -is [uint32] -or
        $Value -is [int64] -or
        $Value -is [uint64]
    ) {
        return $Value.ToString($culture)
    }

    return $null
}

function New-ValueRecord {
    param(
        [AllowNull()]
        [object]$Value,

        [Parameter(Mandatory = $true)]
        [AllowEmptyString()]
        [string]$Unit
    )

    if ([string]::IsNullOrWhiteSpace($Unit)) {
        $Unit = "dimensionless"
    }

    if ($null -eq $Value) {
        return @{
            value      = $null
            unit       = $Unit
            value_text = $null
            value_type = "null"
        }
    }

    $numericText = Get-NumericText $Value

    if (
        ($Value -is [double] -or $Value -is [float]) -and
        (
            [double]::IsNaN([double]$Value) -or
            [double]::IsInfinity([double]$Value)
        )
    ) {
        return @{
            value      = $null
            unit       = $Unit
            value_text = $numericText
            value_type = "non_finite"
        }
    }

    if ($null -ne $numericText) {
        return @{
            value      = [double]$Value
            unit       = $Unit
            value_text = $numericText
            value_type = "number"
        }
    }

    if ($Value -is [bool]) {
        return @{
            value      = [bool]$Value
            unit       = $Unit
            value_text = $Value.ToString()
            value_type = "boolean"
        }
    }

    if ($Value -is [System.Array]) {
        $items = @(
            $Value | ForEach-Object {
                New-ValueRecord $_ "dimensionless"
            }
        )

        return @{
            value      = $items
            unit       = $Unit
            value_text = $null
            value_type = "array"
        }
    }

    return @{
        value      = [string]$Value
        unit       = $Unit
        value_text = [string]$Value
        value_type = "string"
    }
}

function Get-PropertyRecord {
    param(
        [Parameter(Mandatory = $true)]
        [object]$Object,

        [Parameter(Mandatory = $true)]
        [string]$PropertyName
    )

    $unit = "dimensionless"
    $readError = $null

    try {
        $unit = [string][DwsimCaptureReflection]::Invoke(
            $Object,
            "GetPropertyUnit",
            [object[]]@($PropertyName, $null)
        )
    }
    catch {
        $unit = "dimensionless"
    }

    try {
        $value = [DwsimCaptureReflection]::Invoke(
            $Object,
            "GetPropertyValue",
            [object[]]@($PropertyName, $null)
        )

        $record = New-ValueRecord $value $unit
    }
    catch {
        $record = New-ValueRecord $null $unit
        $readError = $_.Exception.Message
    }

    return @{
        property   = $PropertyName
        value      = $record
        unit       = $record.unit
        read_error = $readError
    }
}

function Get-UtilityStates {
    param(
        [AllowNull()]
        [object]$Object
    )

    $states = @()
    $utilities = [DwsimCaptureReflection]::Get(
        $Object,
        "AttachedUtilities"
    )

    foreach ($utility in $utilities) {
        $data = $null
        $readError = $null

        try {
            $data = [DwsimCaptureReflection]::Invoke(
                $utility,
                "SaveData",
                [object[]]@()
            )
        }
        catch {
            $readError = $_.Exception.Message
        }

        $states += @{
            name       = [string][DwsimCaptureReflection]::Get(
                $utility,
                "Name"
            )
            type       = [DwsimCaptureReflection]::TypeName(
                $utility
            )
            data       = $data
            read_error = $readError
        }
    }

    return @(
        $states | Sort-Object name, type
    )
}

function Get-SavedUtilityStates {
    param(
        [Parameter(Mandatory = $true)]
        [string]$CasePath
    )

    if ([IO.Path]::GetExtension($CasePath).Equals(
        ".dwxml",
        [StringComparison]::OrdinalIgnoreCase
    )) {
        [xml]$document = [IO.File]::ReadAllText($CasePath)
    }
    else {
        Add-Type -AssemblyName System.IO.Compression.FileSystem
        $archive = [System.IO.Compression.ZipFile]::OpenRead($CasePath)

        try {
            $entry = @($archive.Entries | Where-Object {
                $_.FullName.EndsWith(".xml")
            }) | Select-Object -First 1

            if ($null -eq $entry) {
                throw "DWSIM case archive contains no XML document"
            }

            $reader = New-Object IO.StreamReader($entry.Open())
            try {
                [xml]$document = $reader.ReadToEnd()
            }
            finally {
                $reader.Dispose()
            }
        }
        finally {
            $archive.Dispose()
        }
    }

    $states = @{}

    foreach ($simulationObject in $document.DWSIM_Simulation_Data.SimulationObjects.SimulationObject) {
        $utilities = @()
        $attachedUtilities = $simulationObject.PSObject.Properties["AttachedUtilities"]
        $utilityEntries = $null

        if ($null -ne $attachedUtilities) {
            $utilityEntries = $attachedUtilities.Value.PSObject.Properties["AttachedUtility"]
        }

        if ($null -ne $utilityEntries) {
            foreach ($utility in $utilityEntries.Value) {
                $data = $null
                $readError = $null

                try {
                    $data = [string]$utility.Data | ConvertFrom-Json
                }
                catch {
                    $readError = $_.Exception.Message
                }

                $utilities += @{
                    name       = [string]$utility.Name
                    type       = [string]$utility.UtilityType
                    data       = $data
                    read_error = $readError
                }
            }
        }

        $states[[string]$simulationObject.Name] = @(
            $utilities | Sort-Object name, type
        )
    }

    return $states
}

function Get-ObjectStates {
    param(
        [Parameter(Mandatory = $true)]
        [object]$Flowsheet,

        [Parameter(Mandatory = $true)]
        [hashtable]$SavedUtilityStates,

        [string[]]$ObjectTags = @(),

        [string[]]$PropertyNames = @(),

        [switch]$CaptureColumnProfiles
    )

    $states = @()

    $simulationObjects = [DwsimCaptureReflection]::Get(
        $Flowsheet,
        "SimulationObjects"
    )

    $objectValues = [DwsimCaptureReflection]::Get(
        $simulationObjects,
        "Values"
    )

    foreach ($object in $objectValues) {
        $name = [DwsimCaptureReflection]::Get(
            $object,
            "Name"
        )

        $graphicObject = [DwsimCaptureReflection]::Get(
            $object,
            "GraphicObject"
        )

        $graphicTag = [DwsimCaptureReflection]::Get(
            $graphicObject,
            "Tag"
        )

        $tag = [string]$name

        if (-not [string]::IsNullOrWhiteSpace(
            [string]$graphicTag
        )) {
            $tag = [string]$graphicTag
        }

        if ($ObjectTags.Count -gt 0 -and $ObjectTags -notcontains $tag) {
            continue
        }

        $properties = @()

        $availablePropertyNames = [DwsimCaptureReflection]::Invoke(
            $object,
            "GetProperties",
            [object[]]@(
                [DWSIM.Interfaces.Enums.PropertyType]::ALL
            )
        )

        foreach ($propertyName in $availablePropertyNames) {
            if ($PropertyNames.Count -gt 0 -and $PropertyNames -notcontains $propertyName) {
                continue
            }

            $properties += Get-PropertyRecord `
                -Object $object `
                -PropertyName ([string]$propertyName)
        }

        $columnProfile = [ordered]@{}

        if ($CaptureColumnProfiles) {
            foreach ($profileName in @("Tf", "Lf", "Vf", "xf", "yf", "Kf")) {
                $profileValue = [DwsimCaptureReflection]::NumericArray(
                    [DwsimCaptureReflection]::Get($object, $profileName)
                )

                if ($null -ne $profileValue) {
                    $columnProfile[$profileName] = $profileValue
                }
            }

            foreach ($profile in @(
                @("Hlf", "xf", "Liquid"),
                @("Hvf", "yf", "Vapor")
            )) {
                $enthalpyProfile = [DwsimCaptureReflection]::ColumnEnthalpyArray(
                    $object,
                    $profile[1],
                    $profile[2]
                )

                if ($null -ne $enthalpyProfile) {
                    $columnProfile[$profile[0]] = $enthalpyProfile
                }
            }

            foreach ($dutyName in @("CondenserDuty", "ReboilerDuty")) {
                $dutyValue = [DwsimCaptureReflection]::Get($object, $dutyName)

                if ($null -ne $dutyValue) {
                    $columnProfile[$dutyName] = New-ValueRecord $dutyValue "kW"
                }
            }
        }

        $states += @{
            tag        = $tag
            name       = [string]$name
            type       = [DwsimCaptureReflection]::TypeName(
                $object
            )
            calculated = [bool][DwsimCaptureReflection]::Get(
                $object,
                "Calculated"
            )
            error      = [DwsimCaptureReflection]::Get(
                $object,
                "ErrorMessage"
            )
            utilities  = @(
                if ($SavedUtilityStates.ContainsKey([string]$name)) {
                    $SavedUtilityStates[[string]$name]
                }
                else {
                    Get-UtilityStates -Object $object
                }
            )
            properties = @(
                $properties | Sort-Object property
            )
            column_profile = $columnProfile
        }
    }

    return @(
        $states | Sort-Object tag
    )
}

function New-CompoundRecord {
    param(
        [Parameter(Mandatory = $true)]
        [string]$CompoundId,

        [Parameter(Mandatory = $true)]
        [object]$Constant,

        [Parameter(Mandatory = $true)]
        [string]$Revision
    )

    $database = Get-MemberValue `
        -Object $Constant `
        -Name "CurrentDB"

    if ([string]::IsNullOrWhiteSpace(
        [string]$database
    )) {
        $database = "DWSIM.AvailableCompounds"
    }

    $canonicalName = [string](
        Get-MemberValue `
            -Object $Constant `
            -Name "Name"
    )

    if ([string]::IsNullOrWhiteSpace($canonicalName)) {
        throw "DWSIM returned an empty canonical name for '$CompoundId'"
    }

    $normalBoilingPoint = [double](
        Get-MemberValue `
            -Object $Constant `
            -Name "NBP"
    )

    $idealHeatCapacity = [DwsimCaptureReflection]::Invoke(
        $Constant,
        "GetIdealGasHeatCapacity",
        [object[]]@(300.0, $null)
    )

    $vaporPressure = [DwsimCaptureReflection]::Invoke(
        $Constant,
        "GetVaporPressure",
        [object[]]@($normalBoilingPoint, $null)
    )

    $elements = [DwsimCaptureReflection]::StringDoubleDictionary(
        (Get-MemberValue -Object $Constant -Name "Elements")
    )

    return @{
        id = $canonicalName

        name = $canonicalName

        cas = [string](
            Get-MemberValue `
                -Object $Constant `
                -Name "CAS_Number"
        )

        formula = [string](
            Get-MemberValue `
                -Object $Constant `
                -Name "Formula"
        )

        elements = $elements

        molecular_weight = New-ValueRecord `
            -Value (
                Get-MemberValue `
                    -Object $Constant `
                    -Name "Molar_Weight"
            ) `
            -Unit "kg/kmol"

        critical_temperature = New-ValueRecord `
            -Value (
                Get-MemberValue `
                    -Object $Constant `
                    -Name "Critical_Temperature"
            ) `
            -Unit "K"

        critical_pressure = New-ValueRecord `
            -Value (
                Get-MemberValue `
                    -Object $Constant `
                    -Name "Critical_Pressure"
            ) `
            -Unit "Pa"

        acentric_factor = New-ValueRecord `
            -Value (
                Get-MemberValue `
                    -Object $Constant `
                    -Name "Acentric_Factor"
            ) `
            -Unit "dimensionless"

        normal_boiling_point = New-ValueRecord `
            -Value $normalBoilingPoint `
            -Unit "K"

        ideal_gas_formation = @{
            temperature = New-ValueRecord `
                -Value 298.15 `
                -Unit "K"

            enthalpy = New-ValueRecord `
                -Value (
                    Get-MemberValue `
                        -Object $Constant `
                        -Name "IG_Enthalpy_of_Formation_25C"
                ) `
                -Unit "kJ/kg"

            gibbs_energy = New-ValueRecord `
                -Value (
                    Get-MemberValue `
                        -Object $Constant `
                        -Name "IG_Gibbs_Energy_of_Formation_25C"
                ) `
                -Unit "kJ/kg"

            entropy = New-ValueRecord `
                -Value (
                    Get-MemberValue `
                        -Object $Constant `
                        -Name "IG_Entropy_of_Formation_25C"
                ) `
                -Unit "kJ/kg/K"
        }

        ideal_reference = @{
            heat_capacity_temperature = New-ValueRecord `
                -Value 300.0 `
                -Unit "K"

            heat_capacity = New-ValueRecord `
                -Value $idealHeatCapacity `
                -Unit "kJ/kg/K"

            vapor_pressure_temperature = New-ValueRecord `
                -Value $normalBoilingPoint `
                -Unit "K"

            vapor_pressure = New-ValueRecord `
                -Value $vaporPressure `
                -Unit "Pa"
        }

        provenance = @{
            database        = [string]$database
            source          = "DWSIM.AvailableCompounds"
            source_revision = $Revision
        }
    }
}

$engineDirectory = (Resolve-Path $EngineBin).Path

$automationPath = Join-Path `
    $engineDirectory `
    "DWSIM.Automation.dll"

$interfacesPath = Join-Path `
    $engineDirectory `
    "DWSIM.Interfaces.dll"

if (-not (Test-Path -LiteralPath $automationPath)) {
    throw "Missing DWSIM.Automation.dll in $engineDirectory"
}

if (-not (Test-Path -LiteralPath $interfacesPath)) {
    throw "Missing DWSIM.Interfaces.dll in $engineDirectory"
}

if ([string]::IsNullOrWhiteSpace($DwsimRevision)) {
    throw "DwsimRevision is required for an auditable capture"
}

$thermoCAssemblyPath = Get-ChildItem `
    -Path $engineDirectory `
    -Filter "ThermoCS.dll" `
    -File `
    -Recurse `
    -ErrorAction SilentlyContinue |
    Select-Object -First 1

[Reflection.Assembly]::LoadFrom(
    $interfacesPath
) | Out-Null

if ($null -ne $thermoCAssemblyPath) {
    [Reflection.Assembly]::LoadFrom(
        $thermoCAssemblyPath.FullName
    ) | Out-Null
}
else {
    Write-Warning (
        "ThermoCS.dll was not found under " +
        "'$engineDirectory'. ThermoC property-package " +
        "discovery may report an error."
    )
}

[Reflection.Assembly]::LoadFrom(
    $automationPath
) | Out-Null

$automation = New-Object `
    DWSIM.Automation.Automation3

$compoundRecords = @()

foreach ($requestedName in $CompoundNames) {
    if (
        -not $automation.AvailableCompounds.ContainsKey(
            $requestedName
        )
    ) {
        throw (
            "Compound '$requestedName' was not found " +
            "in DWSIM.AvailableCompounds"
        )
    }

    $constant = $automation.AvailableCompounds[
        $requestedName
    ]

    $compoundRecords += New-CompoundRecord `
        -CompoundId $requestedName `
        -Constant $constant `
        -Revision $DwsimRevision
}

$caseKind = "compound_catalog"

$caseName = if (
    [string]::IsNullOrWhiteSpace($CaseId)
) {
    "compound-catalog"
}
else {
    $CaseId
}

$inputFile = $null
$inputHash = $null
$before = @()
$after = @()
$propertyPackagesUpdated = 0
$savedUtilityStates = @{}

$solve = @{
    executed = $false
    success  = $true
    errors   = @()
}

if (-not [string]::IsNullOrWhiteSpace($CasePath)) {
    $caseKind = "flowsheet"
    $resolvedCase = (Resolve-Path $CasePath).Path

    $caseName = if (
        [string]::IsNullOrWhiteSpace($CaseId)
    ) {
        [IO.Path]::GetFileNameWithoutExtension(
            $resolvedCase
        )
    }
    else {
        $CaseId
    }

    $inputFile = $resolvedCase

    $inputHash = (
        Get-FileHash `
            -Algorithm SHA256 `
            -Path $resolvedCase
    ).Hash.ToLowerInvariant()

    $savedUtilityStates = Get-SavedUtilityStates `
        -CasePath $resolvedCase

    $flowsheet = $automation.LoadFlowsheet2(
        $resolvedCase
    )

    if ($CalculateBubbleAndDewPoints) {
        $propertyPackagesUpdated = [DwsimCaptureReflection]::SetFlashSetting(
            $flowsheet,
            "CalculateBubbleAndDewPoints",
            "True"
        )
    }

    $before = @(
        Get-ObjectStates `
            -Flowsheet $flowsheet `
            -SavedUtilityStates $savedUtilityStates `
            -ObjectTags $ObjectTags `
            -PropertyNames $PropertyNames `
            -CaptureColumnProfiles:$CaptureColumnProfiles
    )

    $errors = @()

    try {
        $errors = @(
            $automation.CalculateFlowsheet4(
                $flowsheet
            )
        )
    }
    catch {
        $errors = @(
            $_.Exception
        )
    }

    $solve = @{
        executed = $true
        success  = ($errors.Count -eq 0)
        errors   = @(
            $errors | ForEach-Object {
                [string]$_.Message
            }
        )
    }

    $after = @(
        Get-ObjectStates `
            -Flowsheet $flowsheet `
            -SavedUtilityStates $savedUtilityStates `
            -ObjectTags $ObjectTags `
            -PropertyNames $PropertyNames `
            -CaptureColumnProfiles:$CaptureColumnProfiles
    )
}

$document = [ordered]@{
    schema_version = "golden-case-1"
    case_id         = $caseName
    case_kind       = $caseKind

    source = [ordered]@{
        dwsim_revision = $DwsimRevision

        automation_version = [string](
            $automation.GetVersion()
        )

        captured_utc = [DateTime]::UtcNow.ToString(
            "o"
        )

        platform = (
            [Environment]::OSVersion.VersionString
        )

        architecture = (
            [Environment]::GetEnvironmentVariable(
                "PROCESSOR_ARCHITECTURE"
            )
        )

        engine_bin        = $engineDirectory
        input_file        = $inputFile
        input_file_sha256 = $inputHash
        property_package  = $PropertyPackage
        flash_algorithm   = $FlashAlgorithm
        bubble_dew_calculation = [bool]$CalculateBubbleAndDewPoints
        column_profiles_captured = [bool]$CaptureColumnProfiles
        property_packages_updated = $propertyPackagesUpdated

        notes = (
            "Values are captured before and after " +
            "CalculateFlowsheet4; numeric_text preserves " +
            "invariant-culture source formatting."
        )
    }

    inputs = [ordered]@{
        compounds = @(
            $compoundRecords | Sort-Object id
        )

        objects_before = $before
    }

    outputs = [ordered]@{
        solve         = $solve
        objects_after = $after
    }
}

$outputDirectory = Split-Path `
    -Parent `
    $OutputPath

if (
    -not [string]::IsNullOrWhiteSpace(
        $outputDirectory
    ) -and
    -not (Test-Path -LiteralPath $outputDirectory)
) {
    New-Item `
        -ItemType Directory `
        -Path $outputDirectory `
        -Force |
        Out-Null
}

$document |
    ConvertTo-Json -Depth 30 |
    Set-Content `
        -Path $OutputPath `
        -Encoding UTF8

Write-Output "Wrote $OutputPath"
