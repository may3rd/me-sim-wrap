[CmdletBinding()]
param(
    [Parameter(Mandatory=$true)][string]$EngineBin,
    [Parameter(Mandatory=$true)][string]$CasePath,
    [Parameter(Mandatory=$true)][string]$SourceRoot,
    [Parameter(Mandatory=$true)][string]$OutputPath,
    [string]$DwsimRevision="9.0.5.0",
    [double]$TemperatureK=350,
    [double]$PressurePa=1000000,
    [string]$CompositionCsv="0.6,0.4"
)
$ErrorActionPreference="Stop";Set-StrictMode -Version Latest
$composition=@($CompositionCsv.Split(",")|ForEach-Object{[double]::Parse($_,[Globalization.CultureInfo]::InvariantCulture)})
Add-Type -TypeDefinition @"
using System;using System.Collections;using System.Collections.Generic;using System.Globalization;using System.Linq;using System.Reflection;
public static class DwsimBlackOilExtractor{
 static BindingFlags F=BindingFlags.Public|BindingFlags.Instance;
 static object G(object o,string n){var p=o.GetType().GetProperties(F).FirstOrDefault(x=>x.Name==n&&x.GetIndexParameters().Length==0);if(p!=null)return p.GetValue(o,null);var f=o.GetType().GetField(n,F);return f==null?null:f.GetValue(o);}
 static object I(object o,string n,object[] a,Func<MethodInfo,bool> match){return o.GetType().GetMethods(F).First(x=>x.Name==n&&match(x)).Invoke(o,a);}
 static IEnumerable<object> V(object d){foreach(object e in (IEnumerable)d){object v=G(e,"Value");if(v!=null)yield return v;}}
 public static object Package(object fs,string name){foreach(object p in V(G(fs,"PropertyPackages")))if(String.Equals(Convert.ToString(G(p,"ComponentName")),name,StringComparison.Ordinal))return p;throw new InvalidOperationException("Black Oil package not found");}
 public static object Stream(object fs){foreach(object o in V(G(fs,"SimulationObjects")))if(o.GetType().FullName=="DWSIM.Thermodynamics.Streams.MaterialStream")return o;throw new InvalidOperationException("Material stream not found");}
 public static void Current(object pp,object ms){pp.GetType().GetProperties(F).First(x=>x.Name=="CurrentMaterialStream"&&x.CanWrite).SetValue(pp,ms,null);}
 public static object[] Constants(object pp){return ((IEnumerable)I(pp,"DW_GetConstantProperties",new object[0],m=>m.GetParameters().Length==0)).Cast<object>().ToArray();}
 public static string Name(object c){return Convert.ToString(G(c,"Name"),CultureInfo.InvariantCulture);}
 public static double Number(object c,string n){return Convert.ToDouble(G(c,n),CultureInfo.InvariantCulture);}
 public static double VaporPressure(object pp,string name,double t){return Convert.ToDouble(I(pp,"AUX_PVAPi",new object[]{name,t},m=>m.GetParameters().Length==2&&m.GetParameters()[0].ParameterType==typeof(string)),CultureInfo.InvariantCulture);}
 public static double[] VaporizedFractions(object pp,double t,double p){return ((IEnumerable)I(pp,"DW_CalcXY",new object[]{t,p},m=>m.GetParameters().Length==2)).Cast<object>().Select(x=>Convert.ToDouble(x,CultureInfo.InvariantCulture)).ToArray();}
}
"@
$engine=(Resolve-Path -LiteralPath $EngineBin).Path;$case=(Resolve-Path -LiteralPath $CasePath).Path;$source=(Resolve-Path -LiteralPath $SourceRoot).Path
$assembly=Join-Path $engine "DWSIM.Thermodynamics.dll";$packageSource=Join-Path $source "DWSIM.Thermodynamics\PropertyPackages\BlackOil.vb";$modelSource=Join-Path $source "DWSIM.Thermodynamics\PropertyPackages\Models\BlackOilProperties.vb";$flashSource=Join-Path $source "DWSIM.Thermodynamics\FlashAlgorithms\BlackOil.vb";$petroleumSource=Join-Path $source "DWSIM.Thermodynamics\BaseClasses\PropertyMethods.vb"
foreach($f in @((Join-Path $engine "DWSIM.Interfaces.dll"),(Join-Path $engine "DWSIM.Automation.dll"),$assembly,$case,$packageSource,$modelSource,$flashSource,$petroleumSource)){if(-not(Test-Path -LiteralPath $f)){throw "Missing required file: $f"}}
[Reflection.Assembly]::LoadFrom((Join-Path $engine "DWSIM.Interfaces.dll"))|Out-Null;$tc=Get-ChildItem $engine -Filter ThermoCS.dll -File -Recurse -ErrorAction SilentlyContinue|Select-Object -First 1;if($null-ne$tc){[Reflection.Assembly]::LoadFrom($tc.FullName)|Out-Null};[Reflection.Assembly]::LoadFrom($assembly)|Out-Null;[Reflection.Assembly]::LoadFrom((Join-Path $engine "DWSIM.Automation.dll"))|Out-Null
$auto=New-Object DWSIM.Automation.Automation3;$fs=$auto.LoadFlowsheet2($case);$pp=[DwsimBlackOilExtractor]::Package($fs,"Black Oil");[DwsimBlackOilExtractor]::Current($pp,[DwsimBlackOilExtractor]::Stream($fs));$constants=[DwsimBlackOilExtractor]::Constants($pp)
if($constants.Count-ne$composition.Count){throw "Composition length does not match Black Oil case compounds"}
$records=@();$pressures=@();foreach($constant in $constants){$name=[DwsimBlackOilExtractor]::Name($constant);$records+=,[ordered]@{compound_id=$name;specific_gravity_oil=[DwsimBlackOilExtractor]::Number($constant,"BO_SGO");specific_gravity_gas=[DwsimBlackOilExtractor]::Number($constant,"BO_SGG");basic_sediment_water_percent=[DwsimBlackOilExtractor]::Number($constant,"BO_BSW");gas_oil_ratio=[DwsimBlackOilExtractor]::Number($constant,"BO_GOR")};$pressures+=,[DwsimBlackOilExtractor]::VaporPressure($pp,$name,$TemperatureK)}
$vaporized=[DwsimBlackOilExtractor]::VaporizedFractions($pp,$TemperatureK,$PressurePa)
$doc=[ordered]@{schema_version="dwsim-black-oil-data-1";model="Black Oil";source=[ordered]@{product="DWSIM";revision=$DwsimRevision;case_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $case).Hash.ToLowerInvariant();runtime_assembly_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $assembly).Hash.ToLowerInvariant();property_package_source_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $packageSource).Hash.ToLowerInvariant();model_source_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $modelSource).Hash.ToLowerInvariant();flash_source_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $flashSource).Hash.ToLowerInvariant();petroleum_methods_source_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $petroleumSource).Hash.ToLowerInvariant()};compounds=$records;scoped_probe=[ordered]@{temperature_k=$TemperatureK;pressure_pa=$PressurePa;composition=$composition;component_vapor_pressures_pa=$pressures;component_vaporized_fractions=$vaporized}}
$out=[IO.Path]::GetFullPath((Join-Path (Get-Location) $OutputPath));$dir=Split-Path -Parent $out;if($dir){New-Item -ItemType Directory -Force -Path $dir|Out-Null};$doc|ConvertTo-Json -Depth 8|Set-Content -LiteralPath $out -Encoding utf8
[ordered]@{output=$out;component_vapor_pressures_pa=$pressures;component_vaporized_fractions=$vaporized;sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $out).Hash.ToLowerInvariant()}|ConvertTo-Json -Depth 4
