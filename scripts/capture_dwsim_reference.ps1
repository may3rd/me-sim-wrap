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

        PropertyInfo property =
            target.GetType().GetProperty(name, PublicInstance);

        return property == null
            ? null
            : property.GetValue(target, null);
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

function Get-ObjectStates {
    param(
        [Parameter(Mandatory = $true)]
        [object]$Flowsheet
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

        $properties = @()

        $propertyNames = [DwsimCaptureReflection]::Invoke(
            $object,
            "GetProperties",
            [object[]]@(
                [DWSIM.Interfaces.Enums.PropertyType]::ALL
            )
        )

        foreach ($propertyName in $propertyNames) {
            $properties += Get-PropertyRecord `
                -Object $object `
                -PropertyName ([string]$propertyName)
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
            properties = @(
                $properties | Sort-Object property
            )
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
        Get-ObjectStates -Flowsheet $flowsheet
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
        Get-ObjectStates -Flowsheet $flowsheet
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
