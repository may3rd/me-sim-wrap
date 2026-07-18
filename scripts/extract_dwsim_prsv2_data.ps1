[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$DwsimSourceRoot,

    [Parameter(Mandatory = $true)]
    [string]$OutputPath,

    [string]$DwsimRevision = "9.0.5.0"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Convert-InvariantDouble {
    param([Parameter(Mandatory = $true)][string]$Text)
    return [double]::Parse($Text, [Globalization.CultureInfo]::InvariantCulture)
}

$assetDirectory = Join-Path $DwsimSourceRoot "DWSIM.Thermodynamics\Assets"
$alphaPath = Join-Path $assetDirectory "prsv2.dat"
$margulesPath = Join-Path $assetDirectory "prsv2_ip.dat"
$vanLaarPath = Join-Path $assetDirectory "prsv2_ip_vl.dat"
foreach ($path in @($alphaPath, $margulesPath, $vanLaarPath)) {
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Missing DWSIM PRSV2 source asset: $path"
    }
}

$alphaRecords = @(
    Import-Csv -LiteralPath $alphaPath -Delimiter "`t" | ForEach-Object {
        [ordered]@{
            compound = [string]$_.CompoundName
            kappa1 = Convert-InvariantDouble $_.Kappa1
            kappa2 = Convert-InvariantDouble $_.Kappa2
            kappa3 = Convert-InvariantDouble $_.Kappa3
        }
    }
)

function Convert-InteractionTable {
    param([Parameter(Mandatory = $true)][string]$Path)
    return @(
        Import-Csv -LiteralPath $Path -Delimiter "`t" | ForEach-Object {
            [ordered]@{
                first = [string]$_.Compound1
                second = [string]$_.Compound2
                k12 = Convert-InvariantDouble $_.k12
                k21 = Convert-InvariantDouble $_.k21
                reference_temperature_k = Convert-InvariantDouble $_.T_K
            }
        }
    )
}

$margulesRecords = @(Convert-InteractionTable $margulesPath)
$vanLaarRecords = @(Convert-InteractionTable $vanLaarPath)
if ($alphaRecords.Count -ne 90) {
    throw "Expected 90 PRSV2 alpha records, found $($alphaRecords.Count)"
}
if ($margulesRecords.Count -ne 8 -or $vanLaarRecords.Count -ne 8) {
    throw "Expected eight directed PRSV2 interaction pairs per mixing rule"
}

$document = [ordered]@{
    schema_version = "dwsim-prsv2-data-1"
    source = [ordered]@{
        product = "DWSIM"
        revision = $DwsimRevision
        alpha_file = "DWSIM.Thermodynamics/Assets/prsv2.dat"
        alpha_sha256 = (Get-FileHash -LiteralPath $alphaPath -Algorithm SHA256).Hash.ToLowerInvariant()
        margules_file = "DWSIM.Thermodynamics/Assets/prsv2_ip.dat"
        margules_sha256 = (Get-FileHash -LiteralPath $margulesPath -Algorithm SHA256).Hash.ToLowerInvariant()
        van_laar_file = "DWSIM.Thermodynamics/Assets/prsv2_ip_vl.dat"
        van_laar_sha256 = (Get-FileHash -LiteralPath $vanLaarPath -Algorithm SHA256).Hash.ToLowerInvariant()
    }
    alpha_parameters = $alphaRecords
    interactions = [ordered]@{
        margules = $margulesRecords
        van_laar = $vanLaarRecords
    }
}

$resolvedOutput = [IO.Path]::GetFullPath((Join-Path (Get-Location) $OutputPath))
$outputDirectory = Split-Path -Parent $resolvedOutput
if (-not (Test-Path -LiteralPath $outputDirectory -PathType Container)) {
    New-Item -ItemType Directory -Path $outputDirectory -Force | Out-Null
}
$document | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $resolvedOutput -Encoding utf8

[ordered]@{
    output = $resolvedOutput
    alpha_records = $alphaRecords.Count
    margules_pairs = $margulesRecords.Count
    van_laar_pairs = $vanLaarRecords.Count
    sha256 = (Get-FileHash -LiteralPath $resolvedOutput -Algorithm SHA256).Hash.ToLowerInvariant()
} | ConvertTo-Json
