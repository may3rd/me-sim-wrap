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
    [string[]]$CompoundNames = @("Methane", "Ethane", "Propane", "n-Butane", "n-Pentane")
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Get-ClrBaseObject {
    param([AllowNull()][object]$Object)

    if ($null -eq $Object) { return $null }
    return [Management.Automation.PSObject]::AsPSObject($Object).BaseObject
}

function Get-ClrType {
    param([AllowNull()][object]$Object)

    $baseObject = Get-ClrBaseObject $Object
    if ($null -eq $baseObject) { return $null }
    return [object].GetMethod("GetType").Invoke($baseObject, $null)
}

function Get-MemberValue {
    param(
        [AllowNull()]
        [object]$Object,
        [Parameter(Mandatory = $true)]
        [string]$Name
    )

    $baseObject = Get-ClrBaseObject $Object
    if ($null -eq $baseObject) { return $null }
    $flags = [Reflection.BindingFlags]::GetProperty -bor [Reflection.BindingFlags]::Public -bor [Reflection.BindingFlags]::Instance
    try { return (Get-ClrType $baseObject).InvokeMember($Name, $flags, $null, $baseObject, $null) }
    catch [MissingMemberException] { return $null }
}

function Invoke-ClrMethod {
    param(
        [Parameter(Mandatory = $true)][object]$Object,
        [Parameter(Mandatory = $true)][string]$Name,
        [object[]]$Arguments = @()
    )

    $baseObject = Get-ClrBaseObject $Object
    $flags = [Reflection.BindingFlags]::InvokeMethod -bor [Reflection.BindingFlags]::Public -bor [Reflection.BindingFlags]::Instance
    return (Get-ClrType $baseObject).InvokeMember($Name, $flags, $null, $baseObject, $Arguments)
}

function Get-NumericText {
    param([Parameter(Mandatory = $true)][object]$Value)

    $culture = [Globalization.CultureInfo]::InvariantCulture
    if ($Value -is [double]) { return $Value.ToString("R", $culture) }
    if ($Value -is [float]) { return $Value.ToString("R", $culture) }
    if ($Value -is [decimal]) { return $Value.ToString("G29", $culture) }
    if ($Value -is [byte] -or $Value -is [sbyte] -or $Value -is [int16] -or
        $Value -is [uint16] -or $Value -is [int32] -or $Value -is [uint32] -or
        $Value -is [int64] -or $Value -is [uint64]) {
        return $Value.ToString($culture)
    }
    return $null
}

function New-ValueRecord {
    param(
        [AllowNull()]
        [object]$Value,
        [Parameter(Mandatory = $true)]
        [string]$Unit
    )

    if ([string]::IsNullOrWhiteSpace($Unit)) { $Unit = "dimensionless" }

    if ($null -eq $Value) {
        return @{
            value = $null
            unit = $Unit
            value_text = $null
            value_type = "null"
        }
    }

    $numericText = Get-NumericText $Value
    if ($null -ne $numericText) {
        return @{
            value = [double]$Value
            unit = $Unit
            value_text = $numericText
            value_type = "number"
        }
    }

    if ($Value -is [bool]) {
        return @{
            value = [bool]$Value
            unit = $Unit
            value_text = $Value.ToString()
            value_type = "boolean"
        }
    }

    if ($Value -is [System.Array]) {
        $items = @($Value | ForEach-Object { New-ValueRecord $_ "dimensionless" })
        return @{
            value = $items
            unit = $Unit
            value_text = $null
            value_type = "array"
        }
    }

    return @{
        value = [string]$Value
        unit = $Unit
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
    try { $unit = [string](Invoke-ClrMethod -Object $Object -Name "GetPropertyUnit" -Arguments @($PropertyName)) } catch { $unit = "dimensionless" }

    try {
        $value = Invoke-ClrMethod -Object $Object -Name "GetPropertyValue" -Arguments @($PropertyName)
        $record = New-ValueRecord $value $unit
    }
    catch {
        $record = New-ValueRecord $null $unit
        $readError = $_.Exception.Message
    }

    return @{
        property = $PropertyName
        value = $record
        unit = $record.unit
        read_error = $readError
    }
}

function Get-ObjectStates {
    param([Parameter(Mandatory = $true)][object]$Flowsheet)

    $states = @()
    $simulationObjects = Get-MemberValue $Flowsheet "SimulationObjects"
    foreach ($object in (Get-MemberValue $simulationObjects "Values")) {
        $name = Get-MemberValue $object "Name"
        $graphicObject = Get-MemberValue $object "GraphicObject"
        $graphicTag = Get-MemberValue $graphicObject "Tag"
        $tag = [string]$name
        if (-not [string]::IsNullOrWhiteSpace($graphicTag)) {
            $tag = [string]$graphicTag
        }

        $properties = @()
        foreach ($propertyName in (Invoke-ClrMethod -Object $object -Name "GetProperties" -Arguments @([DWSIM.Interfaces.Enums.PropertyType]::ALL))) {
            $properties += Get-PropertyRecord $object ([string]$propertyName)
        }

        $states += @{
            tag = $tag
            name = [string]$name
            type = [string](Get-ClrType $object).Name
            calculated = [bool](Get-MemberValue $object "Calculated")
            error = Get-MemberValue $object "ErrorMessage"
            properties = @($properties | Sort-Object property)
        }
    }
    return @($states | Sort-Object tag)
}

function New-CompoundRecord {
    param(
        [Parameter(Mandatory = $true)]
        [object]$Entry,
        [Parameter(Mandatory = $true)]
        [string]$Revision
    )

    $constant = $Entry.Value
    $database = Get-MemberValue $constant "CurrentDB"
    if ($null -eq $database) { $database = "DWSIM.AvailableCompounds" }

    return @{
        id = [string]$Entry.Key
        name = [string](Get-MemberValue $constant "Name")
        cas = [string](Get-MemberValue $constant "CAS_Number")
        formula = [string](Get-MemberValue $constant "Formula")
        molecular_weight = New-ValueRecord (Get-MemberValue $constant "Molar_Weight") "kg/kmol"
        critical_temperature = New-ValueRecord (Get-MemberValue $constant "Critical_Temperature") "K"
        critical_pressure = New-ValueRecord (Get-MemberValue $constant "Critical_Pressure") "Pa"
        acentric_factor = New-ValueRecord (Get-MemberValue $constant "Acentric_Factor") "dimensionless"
        normal_boiling_point = New-ValueRecord (Get-MemberValue $constant "NBP") "K"
        provenance = @{
            database = [string]$database
            source = "DWSIM.AvailableCompounds"
            source_revision = $Revision
        }
    }
}

$engineDirectory = (Resolve-Path $EngineBin).Path
$automationPath = Join-Path $engineDirectory "DWSIM.Automation.dll"
$interfacesPath = Join-Path $engineDirectory "DWSIM.Interfaces.dll"
$thermoCAssemblyPath = Join-Path $engineDirectory "ThermoCS\ThermoCS.dll"
if (-not (Test-Path $automationPath)) { throw "Missing DWSIM.Automation.dll in $engineDirectory" }
if (-not (Test-Path $interfacesPath)) { throw "Missing DWSIM.Interfaces.dll in $engineDirectory" }
if ([string]::IsNullOrWhiteSpace($DwsimRevision)) { throw "DwsimRevision is required for an auditable capture" }

[Reflection.Assembly]::LoadFrom($interfacesPath) | Out-Null
if (Test-Path $thermoCAssemblyPath) { [Reflection.Assembly]::LoadFrom($thermoCAssemblyPath) | Out-Null }
[Reflection.Assembly]::LoadFrom($automationPath) | Out-Null
$automation = New-Object DWSIM.Automation.Automation3

$compoundRecords = @()
foreach ($requestedName in $CompoundNames) {
    $entry = $automation.AvailableCompounds.GetEnumerator() |
        Where-Object { $_.Key -eq $requestedName } |
        Select-Object -First 1
    if ($null -eq $entry) { throw "Compound '$requestedName' was not found in DWSIM.AvailableCompounds" }
    $compoundRecords += New-CompoundRecord $entry $DwsimRevision
}

$caseKind = "compound_catalog"
$caseName = if ([string]::IsNullOrWhiteSpace($CaseId)) { "compound-catalog" } else { $CaseId }
$inputFile = $null
$inputHash = $null
$before = @()
$after = @()
$solve = @{
    executed = $false
    success = $true
    errors = @()
}

if (-not [string]::IsNullOrWhiteSpace($CasePath)) {
    $caseKind = "flowsheet"
    $resolvedCase = (Resolve-Path $CasePath).Path
    $caseName = if ([string]::IsNullOrWhiteSpace($CaseId)) { [IO.Path]::GetFileNameWithoutExtension($resolvedCase) } else { $CaseId }
    $inputFile = $resolvedCase
    $inputHash = (Get-FileHash -Algorithm SHA256 -Path $resolvedCase).Hash.ToLowerInvariant()
    $flowsheet = $automation.LoadFlowsheet2($resolvedCase)
    $before = @(Get-ObjectStates $flowsheet)

    $errors = @()
    try {
        $errors = @($automation.CalculateFlowsheet4($flowsheet))
    }
    catch {
        $errors = @($_.Exception)
    }

    $solve = @{
        executed = $true
        success = ($errors.Count -eq 0)
        errors = @($errors | ForEach-Object { [string]$_.Message })
    }
    $after = @(Get-ObjectStates $flowsheet)
}

$document = [ordered]@{
    schema_version = "golden-case-1"
    case_id = $caseName
    case_kind = $caseKind
    source = [ordered]@{
        dwsim_revision = $DwsimRevision
        automation_version = [string]$automation.GetVersion()
        captured_utc = [DateTime]::UtcNow.ToString("o")
        platform = [Environment]::OSVersion.VersionString
        architecture = [Environment]::GetEnvironmentVariable("PROCESSOR_ARCHITECTURE")
        engine_bin = $engineDirectory
        input_file = $inputFile
        input_file_sha256 = $inputHash
        property_package = $PropertyPackage
        flash_algorithm = $FlashAlgorithm
        notes = "Values are captured before and after CalculateFlowsheet4; numeric_text preserves invariant-culture source formatting."
    }
    inputs = [ordered]@{
        compounds = @($compoundRecords | Sort-Object id)
        objects_before = $before
    }
    outputs = [ordered]@{
        solve = $solve
        objects_after = $after
    }
}

$outputDirectory = Split-Path -Parent $OutputPath
if (-not [string]::IsNullOrWhiteSpace($outputDirectory) -and -not (Test-Path $outputDirectory)) {
    New-Item -ItemType Directory -Path $outputDirectory -Force | Out-Null
}

$document | ConvertTo-Json -Depth 30 | Set-Content -Path $OutputPath -Encoding UTF8
Write-Output "Wrote $OutputPath"
