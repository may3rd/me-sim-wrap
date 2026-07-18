[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$EngineBin,
    [Parameter(Mandatory = $true)][string]$CasePath,
    [Parameter(Mandatory = $true)][string]$SourceRoot,
    [Parameter(Mandatory = $true)][string]$OutputPath,
    [string]$DwsimRevision = "9.0.5.0"
)
$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

Add-Type -TypeDefinition @"
using System;using System.Collections;using System.Globalization;using System.Linq;using System.Reflection;
public static class DwsimChaoSeaderExtractor {
 static BindingFlags F=BindingFlags.Public|BindingFlags.Instance;
 static object G(object o,string n){var p=o.GetType().GetProperties(F).FirstOrDefault(x=>x.Name==n&&x.GetIndexParameters().Length==0);if(p!=null)return p.GetValue(o,null);var f=o.GetType().GetField(n,F);return f==null?null:f.GetValue(o);}
 static object I(object o,string n){return o.GetType().GetMethods(F).First(x=>x.Name==n&&x.GetParameters().Length==0).Invoke(o,null);}
 public static object Package(object fs){foreach(DictionaryEntry e in(IDictionary)G(fs,"PropertyPackages"))if(String.Equals(Convert.ToString(G(e.Value,"ComponentName"),CultureInfo.InvariantCulture),"Chao-Seader",StringComparison.Ordinal))return e.Value;throw new InvalidOperationException("No exact Chao-Seader package match.");}
 public static void SetStream(object fs,object pp){object ms=null;foreach(DictionaryEntry e in(IDictionary)G(fs,"SimulationObjects"))if(e.Value!=null&&e.Value.GetType().FullName=="DWSIM.Thermodynamics.Streams.MaterialStream")ms=e.Value;if(ms==null)throw new InvalidOperationException("No material stream.");pp.GetType().GetProperties(F).First(x=>x.Name=="CurrentMaterialStream"&&x.CanWrite).SetValue(pp,ms,null);}
 public static string[] S(object pp,string n){return((IEnumerable)I(pp,n)).Cast<object>().Select(x=>Convert.ToString(x,CultureInfo.InvariantCulture)).ToArray();}
 public static double[] D(object pp,string n){return((IEnumerable)I(pp,n)).Cast<object>().Select(x=>Convert.ToDouble(x,CultureInfo.InvariantCulture)).ToArray();}
 public static string TypeName(object pp){return pp.GetType().FullName;}
}
"@

$engine=(Resolve-Path -LiteralPath $EngineBin).Path
$case=(Resolve-Path -LiteralPath $CasePath).Path
$source=(Resolve-Path -LiteralPath $SourceRoot).Path
$assembly=Join-Path $engine "DWSIM.Thermodynamics.dll"
$packageSource=Join-Path $source "DWSIM.Thermodynamics\PropertyPackages\ChaoSeader.vb"
$modelSource=Join-Path $source "DWSIM.Thermodynamics\PropertyPackages\Models\ChaoSeader.vb"
foreach($f in @("DWSIM.Interfaces.dll","DWSIM.Automation.dll","DWSIM.Thermodynamics.dll")){if(-not(Test-Path -LiteralPath (Join-Path $engine $f))){throw "Missing $f"}}
foreach($f in @($packageSource,$modelSource)){if(-not(Test-Path -LiteralPath $f)){throw "Missing source file: $f"}}
[Reflection.Assembly]::LoadFrom((Join-Path $engine "DWSIM.Interfaces.dll"))|Out-Null
$tc=Get-ChildItem $engine -Filter ThermoCS.dll -File -Recurse -ErrorAction SilentlyContinue|Select-Object -First 1;if($null-ne$tc){[Reflection.Assembly]::LoadFrom($tc.FullName)|Out-Null}
[Reflection.Assembly]::LoadFrom($assembly)|Out-Null
[Reflection.Assembly]::LoadFrom((Join-Path $engine "DWSIM.Automation.dll"))|Out-Null
$auto=New-Object DWSIM.Automation.Automation3
$fs=$auto.LoadFlowsheet2($case)
$pp=[DwsimChaoSeaderExtractor]::Package($fs)
[DwsimChaoSeaderExtractor]::SetStream($fs,$pp)
$names=[DwsimChaoSeaderExtractor]::S($pp,"RET_VNAMES")
$mw=[DwsimChaoSeaderExtractor]::D($pp,"RET_VMM")
$tcv=[DwsimChaoSeaderExtractor]::D($pp,"RET_VTC")
$pc=[DwsimChaoSeaderExtractor]::D($pp,"RET_VPC")
$w=[DwsimChaoSeaderExtractor]::D($pp,"RET_VW")
$vl=[DwsimChaoSeaderExtractor]::D($pp,"RET_VVL")
$csa=[DwsimChaoSeaderExtractor]::D($pp,"RET_VCSAc")
$css=[DwsimChaoSeaderExtractor]::D($pp,"RET_VCSS")
$records=@(for($i=0;$i-lt$names.Count;$i++){[ordered]@{compound_id=$names[$i];molecular_weight=$mw[$i];critical_temperature_k=$tcv[$i];critical_pressure_pa=$pc[$i];acentric_factor=$w[$i];liquid_molar_volume=$vl[$i];chao_seader_acentricity=$csa[$i];solubility_parameter=$css[$i]}})
$doc=[ordered]@{schema_version="dwsim-chao-seader-data-1";model="Chao-Seader";source=[ordered]@{product="DWSIM";revision=$DwsimRevision;case_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $case).Hash.ToLowerInvariant();runtime_assembly_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $assembly).Hash.ToLowerInvariant();property_package_class=[DwsimChaoSeaderExtractor]::TypeName($pp);property_package_source="DWSIM.Thermodynamics/PropertyPackages/ChaoSeader.vb";property_package_source_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $packageSource).Hash.ToLowerInvariant();model_source="DWSIM.Thermodynamics/PropertyPackages/Models/ChaoSeader.vb";model_source_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $modelSource).Hash.ToLowerInvariant()};compounds=$records}
$out=[IO.Path]::GetFullPath((Join-Path (Get-Location) $OutputPath));$dir=Split-Path -Parent $out;if($dir){New-Item -ItemType Directory -Force -Path $dir|Out-Null}
$doc|ConvertTo-Json -Depth 8|Set-Content -LiteralPath $out -Encoding utf8
[ordered]@{output=$out;compounds=$records.Count;sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $out).Hash.ToLowerInvariant()}|ConvertTo-Json
