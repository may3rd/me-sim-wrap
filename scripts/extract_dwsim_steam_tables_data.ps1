[CmdletBinding()]
param([Parameter(Mandatory=$true)][string]$EngineBin,[Parameter(Mandatory=$true)][string]$CasePath,[Parameter(Mandatory=$true)][string]$SourceRoot,[Parameter(Mandatory=$true)][string]$OutputPath,[string]$DwsimRevision="9.0.5.0")
$ErrorActionPreference="Stop";Set-StrictMode -Version Latest
$engine=(Resolve-Path -LiteralPath $EngineBin).Path;$case=(Resolve-Path -LiteralPath $CasePath).Path;$source=(Resolve-Path -LiteralPath $SourceRoot).Path
$assembly=Join-Path $engine "DWSIM.Thermodynamics.dll";$packageSource=Join-Path $source "DWSIM.Thermodynamics\PropertyPackages\SteamTables.vb";$modelSource=Join-Path $source "DWSIM.Thermodynamics\PropertyPackages\Models\SteamTables_IAPWS_IF97.vb"
foreach($f in @($assembly,$case,$packageSource,$modelSource)){if(-not(Test-Path -LiteralPath $f)){throw "Missing required file: $f"}}
$text=Get-Content -LiteralPath $modelSource -Raw
$coefficients=@(for($i=1;$i-le 10;$i++){$match=[regex]::Match($text,"nreg4\("+$i+"\)\s*=\s*([-+0-9.Ee]+)");if(-not$match.Success){throw "Missing IAPWS region-4 coefficient $i"};[double]::Parse($match.Groups[1].Value,[Globalization.CultureInfo]::InvariantCulture)})
$doc=[ordered]@{schema_version="dwsim-steam-tables-data-1";model="Steam Tables";source=[ordered]@{product="DWSIM";revision=$DwsimRevision;case_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $case).Hash.ToLowerInvariant();runtime_assembly_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $assembly).Hash.ToLowerInvariant();property_package_source_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $packageSource).Hash.ToLowerInvariant();iapws_source_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $modelSource).Hash.ToLowerInvariant()};temperature_range_k=@(273.15,647.096);region4_coefficients=$coefficients}
$out=[IO.Path]::GetFullPath((Join-Path (Get-Location) $OutputPath));$dir=Split-Path -Parent $out;if($dir){New-Item -ItemType Directory -Force -Path $dir|Out-Null};$doc|ConvertTo-Json -Depth 6|Set-Content -LiteralPath $out -Encoding utf8
[ordered]@{output=$out;coefficients=$coefficients.Count;sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $out).Hash.ToLowerInvariant()}|ConvertTo-Json
