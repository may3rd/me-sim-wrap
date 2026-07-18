[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$EngineBin,
    [Parameter(Mandatory = $true)][string]$CasePath,
    [Parameter(Mandatory = $true)][string]$OutputPath,
    [string]$DwsimRevision = "9.0.5.0",
    [string]$PropertyPackageName = "Modified UNIFAC (Dortmund)",
    [string]$Model = "Modified UNIFAC (Dortmund)",
    [string]$GroupsResource = "DWSIM.Thermodynamics.modfac.txt",
    [string]$InteractionsResource = "DWSIM.Thermodynamics.modfac_ip.txt"
)
$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

Add-Type -TypeDefinition @"
using System;using System.Collections;using System.Collections.Generic;using System.Globalization;
using System.IO;using System.Linq;using System.Reflection;using System.Security.Cryptography;
public sealed class MFGroup { public int P,S; public string PN,N; public double R,Q; }
public sealed class MFPair { public int F,S; public double A12,B12,C12,A21,B21,C21; }
public sealed class MFFraction { public int S; public double V; }
public sealed class MFBasis { public string N; public double Q,R; public MFFraction[] E; }
public static class DwsimModfacExtractor {
 static BindingFlags F=BindingFlags.Public|BindingFlags.Instance;
 static object G(object o,string n){var p=o.GetType().GetProperties(F).FirstOrDefault(x=>x.Name==n&&x.GetIndexParameters().Length==0);if(p!=null)return p.GetValue(o,null);var f=o.GetType().GetField(n,F);return f==null?null:f.GetValue(o);}
 static object I(object o,string n){return o.GetType().GetMethods(F).First(x=>x.Name==n&&x.GetParameters().Length==0).Invoke(o,null);}
 public static object Package(object fs,string name){foreach(DictionaryEntry e in(IDictionary)G(fs,"PropertyPackages"))if(String.Equals(Convert.ToString(G(e.Value,"ComponentName"),CultureInfo.InvariantCulture),name,StringComparison.Ordinal))return e.Value;throw new InvalidOperationException("No exact modified-UNIFAC package match: "+name);}
 public static void SetStream(object fs,object pp){object ms=null;foreach(DictionaryEntry e in(IDictionary)G(fs,"SimulationObjects"))if(e.Value!=null&&e.Value.GetType().FullName=="DWSIM.Thermodynamics.Streams.MaterialStream")ms=e.Value;if(ms==null)throw new InvalidOperationException("No material stream.");pp.GetType().GetProperties(F).First(x=>x.Name=="CurrentMaterialStream"&&x.CanWrite).SetValue(pp,ms,null);}
 static double[] D(object pp,string n){return((IEnumerable)I(pp,n)).Cast<object>().Select(x=>Convert.ToDouble(x,CultureInfo.InvariantCulture)).ToArray();}
 public static MFBasis[] Basis(object pp){var names=((IEnumerable)I(pp,"RET_VNAMES")).Cast<object>().Select(Convert.ToString).ToArray();var q=D(pp,"RET_VQ");var r=D(pp,"RET_VR");var eki=((IEnumerable)I(pp,"RET_VEKI")).Cast<object>().ToArray();return names.Select((n,i)=>{var e=new List<MFFraction>();foreach(DictionaryEntry x in(IDictionary)eki[i])e.Add(new MFFraction{S=Convert.ToInt32(x.Key),V=Convert.ToDouble(x.Value,CultureInfo.InvariantCulture)});return new MFBasis{N=n,Q=q[i],R=r[i],E=e.OrderBy(x=>x.S).ToArray()};}).ToArray();}
 static object GS(object pp){return G(I(pp,"GetModel"),"ModfGroups");}
 public static MFGroup[] Groups(object pp){var rows=new List<MFGroup>();foreach(DictionaryEntry e in(IDictionary)G(GS(pp),"Groups")){var x=e.Value;rows.Add(new MFGroup{P=Convert.ToInt32(G(x,"PrimaryGroup")),S=Convert.ToInt32(G(x,"Secondary_Group")),PN=Convert.ToString(G(x,"MainGroupName")),N=Convert.ToString(G(x,"GroupName")),R=Convert.ToDouble(G(x,"R")),Q=Convert.ToDouble(G(x,"Q"))});}return rows.OrderBy(x=>x.S).ToArray();}
 static double V(IDictionary d,object a,object b){return Convert.ToDouble(((IDictionary)d[a])[b],CultureInfo.InvariantCulture);}
 public static MFPair[] Pairs(object pp){var g=GS(pp);var a12=(IDictionary)G(g,"InteracParam_aij");var b12=(IDictionary)G(g,"InteracParam_bij");var c12=(IDictionary)G(g,"InteracParam_cij");var a21=(IDictionary)G(g,"InteracParam_aji");var b21=(IDictionary)G(g,"InteracParam_bji");var c21=(IDictionary)G(g,"InteracParam_cji");var rows=new List<MFPair>();foreach(DictionaryEntry x in a12)foreach(DictionaryEntry y in(IDictionary)x.Value)rows.Add(new MFPair{F=Convert.ToInt32(x.Key),S=Convert.ToInt32(y.Key),A12=V(a12,x.Key,y.Key),B12=V(b12,x.Key,y.Key),C12=V(c12,x.Key,y.Key),A21=a21==null?0:V(a21,x.Key,y.Key),B21=b21==null?0:V(b21,x.Key,y.Key),C21=c21==null?0:V(c21,x.Key,y.Key)});return rows.OrderBy(x=>x.F).ThenBy(x=>x.S).ToArray();}
 public static string Hash(object pp,string n){using(Stream s=pp.GetType().Assembly.GetManifestResourceStream(n))using(SHA256 h=SHA256.Create())return String.Concat(h.ComputeHash(s).Select(x=>x.ToString("x2")));}
}
"@
$engine=(Resolve-Path -LiteralPath $EngineBin).Path;$case=(Resolve-Path -LiteralPath $CasePath).Path
foreach($f in @("DWSIM.Interfaces.dll","DWSIM.Thermodynamics.dll","DWSIM.Automation.dll")){if(-not(Test-Path (Join-Path $engine $f))){throw "Missing $f"}}
[Reflection.Assembly]::LoadFrom((Join-Path $engine "DWSIM.Interfaces.dll"))|Out-Null
$tc=Get-ChildItem $engine -Filter ThermoCS.dll -File -Recurse -ErrorAction SilentlyContinue|Select-Object -First 1;if($null-ne$tc){[Reflection.Assembly]::LoadFrom($tc.FullName)|Out-Null}
[Reflection.Assembly]::LoadFrom((Join-Path $engine "DWSIM.Thermodynamics.dll"))|Out-Null;[Reflection.Assembly]::LoadFrom((Join-Path $engine "DWSIM.Automation.dll"))|Out-Null
$auto=New-Object DWSIM.Automation.Automation3;$fs=$auto.LoadFlowsheet2($case);$pp=[DwsimModfacExtractor]::Package($fs,$PropertyPackageName);[DwsimModfacExtractor]::SetStream($fs,$pp)
$basis=@([DwsimModfacExtractor]::Basis($pp)|ForEach-Object{[ordered]@{compound_id=$_.N;q=$_.Q;r=$_.R;group_surface_fractions=@($_.E|ForEach-Object{[ordered]@{secondary_id=$_.S;value=$_.V}})}})
$groups=@([DwsimModfacExtractor]::Groups($pp)|ForEach-Object{[ordered]@{primary_id=$_.P;secondary_id=$_.S;primary_name=$_.PN;group_name=$_.N;r=$_.R;q=$_.Q}})
$pairs=@([DwsimModfacExtractor]::Pairs($pp)|ForEach-Object{[ordered]@{first_primary_id=$_.F;second_primary_id=$_.S;a12=$_.A12;b12=$_.B12;c12=$_.C12;a21=$_.A21;b21=$_.B21;c21=$_.C21}})
if($basis.Count-ne2-or$groups.Count-eq0-or$pairs.Count-eq0){throw "Invalid Dortmund extraction count"}
$doc=[ordered]@{schema_version="dwsim-modfac-data-1";model=$Model;source=[ordered]@{product="DWSIM";revision=$DwsimRevision;case_sha256=(Get-FileHash $case -Algorithm SHA256).Hash.ToLowerInvariant();groups_resource=$GroupsResource;groups_sha256=[DwsimModfacExtractor]::Hash($pp,$GroupsResource);interactions_resource=$InteractionsResource;interactions_sha256=[DwsimModfacExtractor]::Hash($pp,$InteractionsResource)};compound_basis=$basis;groups=$groups;interaction_pairs=$pairs}
$out=[IO.Path]::GetFullPath((Join-Path (Get-Location) $OutputPath));$dir=Split-Path -Parent $out;if(-not(Test-Path $dir)){New-Item -ItemType Directory $dir|Out-Null};$doc|ConvertTo-Json -Depth 10|Set-Content $out -Encoding utf8
[ordered]@{output=$out;groups=$groups.Count;interaction_pairs=$pairs.Count;compound_basis_records=$basis.Count;sha256=(Get-FileHash $out -Algorithm SHA256).Hash.ToLowerInvariant()}|ConvertTo-Json
