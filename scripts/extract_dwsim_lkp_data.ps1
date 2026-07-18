[CmdletBinding()]
param([Parameter(Mandatory=$true)][string]$EngineBin,[Parameter(Mandatory=$true)][string]$CasePath,[Parameter(Mandatory=$true)][string]$SourceRoot,[Parameter(Mandatory=$true)][string]$OutputPath,[string]$DwsimRevision="9.0.5.0")
$ErrorActionPreference="Stop";Set-StrictMode -Version Latest
Add-Type -TypeDefinition @"
using System;using System.Collections;using System.Collections.Generic;using System.Globalization;using System.Linq;using System.Reflection;using System.IO;using System.Security.Cryptography;
public sealed class LKPC {public string N;public double M,T,P,W,V;}
public sealed class LKPI {public string A,B;public double K;}
public static class DwsimLKPExtractor {
 static BindingFlags F=BindingFlags.Public|BindingFlags.Instance;
 static object G(object o,string n){var p=o.GetType().GetProperties(F).FirstOrDefault(x=>x.Name==n&&x.GetIndexParameters().Length==0);if(p!=null)return p.GetValue(o,null);var f=o.GetType().GetField(n,F);return f==null?null:f.GetValue(o);}
 static object I(object o,string n){return o.GetType().GetMethods(F).First(x=>x.Name==n&&x.GetParameters().Length==0).Invoke(o,null);}
 static double[] D(object pp,string n){return((IEnumerable)I(pp,n)).Cast<object>().Select(x=>Convert.ToDouble(x,CultureInfo.InvariantCulture)).ToArray();}
 public static object Package(object fs){foreach(DictionaryEntry e in(IDictionary)G(fs,"PropertyPackages"))if(e.Value!=null&&e.Value.GetType().FullName=="DWSIM.Thermodynamics.PropertyPackages.LKPPropertyPackage")return e.Value;throw new InvalidOperationException("No exact LKP package class match.");}
 public static void SetStream(object fs,object pp){object ms=null;foreach(DictionaryEntry e in(IDictionary)G(fs,"SimulationObjects"))if(e.Value!=null&&e.Value.GetType().FullName=="DWSIM.Thermodynamics.Streams.MaterialStream")ms=e.Value;if(ms==null)throw new InvalidOperationException("No material stream.");pp.GetType().GetProperties(F).First(x=>x.Name=="CurrentMaterialStream"&&x.CanWrite).SetValue(pp,ms,null);}
 public static LKPC[] Compounds(object pp){var n=((IEnumerable)I(pp,"RET_VNAMES")).Cast<object>().Select(Convert.ToString).ToArray();var m=D(pp,"RET_VMM");var t=D(pp,"RET_VTC");var p=D(pp,"RET_VPC");var w=D(pp,"RET_VW");var v=D(pp,"RET_VVC");return n.Select((x,i)=>new LKPC{N=x,M=m[i],T=t[i],P=p[i],W=w[i],V=v[i]}).ToArray();}
 public static LKPI[] Pairs(object pp){var model=I(pp,"GetModel");var rows=new List<LKPI>();foreach(DictionaryEntry a in(IDictionary)G(model,"InteractionParameters"))foreach(DictionaryEntry b in(IDictionary)a.Value)rows.Add(new LKPI{A=Convert.ToString(a.Key),B=Convert.ToString(b.Key),K=Convert.ToDouble(G(b.Value,"kij"),CultureInfo.InvariantCulture)});return rows.OrderBy(x=>x.A,StringComparer.Ordinal).ThenBy(x=>x.B,StringComparer.Ordinal).ToArray();}
 public static string Hash(object pp,string name){using(Stream s=pp.GetType().Assembly.GetManifestResourceStream(name))using(SHA256 h=SHA256.Create())return String.Concat(h.ComputeHash(s).Select(x=>x.ToString("x2")));}
 public static string TypeName(object pp){return pp.GetType().FullName;}
 public static string Name(object pp){return Convert.ToString(G(pp,"ComponentName"),CultureInfo.InvariantCulture);}
}
"@
$engine=(Resolve-Path -LiteralPath $EngineBin).Path;$case=(Resolve-Path -LiteralPath $CasePath).Path;$source=(Resolve-Path -LiteralPath $SourceRoot).Path
$assembly=Join-Path $engine "DWSIM.Thermodynamics.dll";$packageSource=Join-Path $source "DWSIM.Thermodynamics\PropertyPackages\LeeKeslerPlocker.vb";$modelSource=Join-Path $source "DWSIM.Thermodynamics\PropertyPackages\Models\LeeKeslerPlocker.vb"
foreach($f in @((Join-Path $engine "DWSIM.Interfaces.dll"),(Join-Path $engine "DWSIM.Automation.dll"),$assembly,$packageSource,$modelSource)){if(-not(Test-Path -LiteralPath $f)){throw "Missing required file: $f"}}
[Reflection.Assembly]::LoadFrom((Join-Path $engine "DWSIM.Interfaces.dll"))|Out-Null;$tc=Get-ChildItem $engine -Filter ThermoCS.dll -File -Recurse -ErrorAction SilentlyContinue|Select-Object -First 1;if($null-ne$tc){[Reflection.Assembly]::LoadFrom($tc.FullName)|Out-Null};[Reflection.Assembly]::LoadFrom($assembly)|Out-Null;[Reflection.Assembly]::LoadFrom((Join-Path $engine "DWSIM.Automation.dll"))|Out-Null
$auto=New-Object DWSIM.Automation.Automation3;$fs=$auto.LoadFlowsheet2($case);$pp=[DwsimLKPExtractor]::Package($fs);[DwsimLKPExtractor]::SetStream($fs,$pp)
$compounds=@([DwsimLKPExtractor]::Compounds($pp)|ForEach-Object{[ordered]@{compound_id=$_.N;molecular_weight=$_.M;critical_temperature_k=$_.T;critical_pressure_pa=$_.P;acentric_factor=$_.W;critical_volume=$_.V}})
$pairs=@([DwsimLKPExtractor]::Pairs($pp)|ForEach-Object{[ordered]@{first_compound_id=$_.A;second_compound_id=$_.B;kij=$_.K}})
$doc=[ordered]@{schema_version="dwsim-lkp-data-1";model=[DwsimLKPExtractor]::Name($pp);source=[ordered]@{product="DWSIM";revision=$DwsimRevision;case_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $case).Hash.ToLowerInvariant();runtime_assembly_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $assembly).Hash.ToLowerInvariant();interaction_resource="DWSIM.Thermodynamics.lkp_ip.dat";interaction_resource_sha256=[DwsimLKPExtractor]::Hash($pp,"DWSIM.Thermodynamics.lkp_ip.dat");property_package_class=[DwsimLKPExtractor]::TypeName($pp);property_package_source_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $packageSource).Hash.ToLowerInvariant();model_source_sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $modelSource).Hash.ToLowerInvariant()};compounds=$compounds;interaction_pairs=$pairs}
$out=[IO.Path]::GetFullPath((Join-Path (Get-Location) $OutputPath));$dir=Split-Path -Parent $out;if($dir){New-Item -ItemType Directory -Force -Path $dir|Out-Null};$doc|ConvertTo-Json -Depth 8|Set-Content -LiteralPath $out -Encoding utf8
[ordered]@{output=$out;compounds=$compounds.Count;interaction_pairs=$pairs.Count;sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $out).Hash.ToLowerInvariant()}|ConvertTo-Json
