[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$EngineBin,
    [Parameter(Mandatory = $true)][string]$CasePath,
    [Parameter(Mandatory = $true)][string]$OutputPath,
    [string]$DwsimRevision = "9.0.5.0"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

Add-Type -TypeDefinition @"
using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Security.Cryptography;

public sealed class UnifacGroupRow
{
    public int PrimaryId { get; set; }
    public int SecondaryId { get; set; }
    public string PrimaryName { get; set; }
    public string GroupName { get; set; }
    public double R { get; set; }
    public double Q { get; set; }
}
public sealed class UnifacInteractionRow
{
    public int FirstPrimaryId { get; set; }
    public int SecondPrimaryId { get; set; }
    public double A { get; set; }
}
public sealed class UnifacSurfaceFraction
{
    public int SecondaryId { get; set; }
    public double Value { get; set; }
}
public sealed class UnifacCompoundBasis
{
    public string Name { get; set; }
    public double Q { get; set; }
    public double R { get; set; }
    public UnifacSurfaceFraction[] SurfaceFractions { get; set; }
}

public static class DwsimUnifacExtractor
{
    private const BindingFlags PublicInstance = BindingFlags.Public | BindingFlags.Instance;
    private static object Get(object target, string name)
    {
        PropertyInfo property = target.GetType().GetProperties(PublicInstance)
            .FirstOrDefault(candidate => candidate.Name == name && candidate.GetIndexParameters().Length == 0);
        if (property != null) return property.GetValue(target, null);
        FieldInfo field = target.GetType().GetField(name, PublicInstance);
        return field == null ? null : field.GetValue(target);
    }
    private static object Invoke(object target, string name)
    {
        MethodInfo method = target.GetType().GetMethods(PublicInstance)
            .First(candidate => candidate.Name == name && candidate.GetParameters().Length == 0);
        return method.Invoke(target, null);
    }
    public static object PropertyPackage(object flowsheet)
    {
        foreach (DictionaryEntry entry in (IDictionary)Get(flowsheet, "PropertyPackages"))
            if (String.Equals(Convert.ToString(Get(entry.Value, "ComponentName"),
                    CultureInfo.InvariantCulture), "UNIFAC", StringComparison.Ordinal)) return entry.Value;
        throw new InvalidOperationException("No exact UNIFAC property-package match was found.");
    }
    public static void SetCurrentMaterialStream(object flowsheet, object package)
    {
        object stream = null;
        foreach (DictionaryEntry entry in (IDictionary)Get(flowsheet, "SimulationObjects"))
            if (entry.Value != null && entry.Value.GetType().FullName ==
                    "DWSIM.Thermodynamics.Streams.MaterialStream") stream = entry.Value;
        if (stream == null) throw new InvalidOperationException("The case has no material stream.");
        package.GetType().GetProperties(PublicInstance)
            .First(candidate => candidate.Name == "CurrentMaterialStream" && candidate.CanWrite)
            .SetValue(package, stream, null);
    }
    private static double[] Values(object package, string method)
    {
        return ((IEnumerable)Invoke(package, method)).Cast<object>()
            .Select(value => Convert.ToDouble(value, CultureInfo.InvariantCulture)).ToArray();
    }
    public static UnifacCompoundBasis[] CompoundBasis(object package)
    {
        string[] names = ((IEnumerable)Invoke(package, "RET_VNAMES")).Cast<object>()
            .Select(value => Convert.ToString(value, CultureInfo.InvariantCulture)).ToArray();
        double[] q = Values(package, "RET_VQ");
        double[] r = Values(package, "RET_VR");
        object[] eki = ((IEnumerable)Invoke(package, "RET_VEKI")).Cast<object>().ToArray();
        return names.Select((name, index) => {
            var fractions = new List<UnifacSurfaceFraction>();
            foreach (DictionaryEntry entry in (IDictionary)eki[index])
                fractions.Add(new UnifacSurfaceFraction {
                    SecondaryId = Convert.ToInt32(entry.Key, CultureInfo.InvariantCulture),
                    Value = Convert.ToDouble(entry.Value, CultureInfo.InvariantCulture)
                });
            return new UnifacCompoundBasis {
                Name = name, Q = q[index], R = r[index],
                SurfaceFractions = fractions.OrderBy(value => value.SecondaryId).ToArray()
            };
        }).ToArray();
    }
    private static object Groups(object package)
    {
        object model = Invoke(package, "GetModel");
        return Get(model, "UnifGroups");
    }
    public static UnifacGroupRow[] GroupRows(object package)
    {
        IDictionary groups = (IDictionary)Get(Groups(package), "Groups");
        var rows = new List<UnifacGroupRow>();
        foreach (DictionaryEntry entry in groups)
        {
            object group = entry.Value;
            rows.Add(new UnifacGroupRow {
                PrimaryId = Convert.ToInt32(Get(group, "PrimaryGroup"), CultureInfo.InvariantCulture),
                SecondaryId = Convert.ToInt32(Get(group, "Secondary_Group"), CultureInfo.InvariantCulture),
                PrimaryName = Convert.ToString(Get(group, "PrimGroupName"), CultureInfo.InvariantCulture),
                GroupName = Convert.ToString(Get(group, "GroupName"), CultureInfo.InvariantCulture),
                R = Convert.ToDouble(Get(group, "R"), CultureInfo.InvariantCulture),
                Q = Convert.ToDouble(Get(group, "Q"), CultureInfo.InvariantCulture)
            });
        }
        return rows.OrderBy(row => row.SecondaryId).ToArray();
    }
    public static UnifacInteractionRow[] InteractionRows(object package)
    {
        IDictionary outer = (IDictionary)Get(Groups(package), "InteracParam");
        var rows = new List<UnifacInteractionRow>();
        foreach (DictionaryEntry first in outer)
            foreach (DictionaryEntry second in (IDictionary)first.Value)
                rows.Add(new UnifacInteractionRow {
                    FirstPrimaryId = Convert.ToInt32(first.Key, CultureInfo.InvariantCulture),
                    SecondPrimaryId = Convert.ToInt32(second.Key, CultureInfo.InvariantCulture),
                    A = Convert.ToDouble(second.Value, CultureInfo.InvariantCulture)
                });
        return rows.OrderBy(row => row.FirstPrimaryId).ThenBy(row => row.SecondPrimaryId).ToArray();
    }
    public static string ResourceSha256(object package, string resourceName)
    {
        using (Stream stream = package.GetType().Assembly.GetManifestResourceStream(resourceName))
        using (SHA256 sha = SHA256.Create())
            return String.Concat(sha.ComputeHash(stream).Select(value => value.ToString("x2")));
    }
}
"@

$engineDirectory = (Resolve-Path -LiteralPath $EngineBin).Path
$resolvedCase = (Resolve-Path -LiteralPath $CasePath).Path
foreach ($file in @("DWSIM.Interfaces.dll", "DWSIM.Thermodynamics.dll", "DWSIM.Automation.dll")) {
    $path = Join-Path $engineDirectory $file
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) { throw "Missing DWSIM assembly: $path" }
}
[Reflection.Assembly]::LoadFrom((Join-Path $engineDirectory "DWSIM.Interfaces.dll")) | Out-Null
$thermoC = Get-ChildItem -Path $engineDirectory -Filter "ThermoCS.dll" -File -Recurse `
    -ErrorAction SilentlyContinue | Select-Object -First 1
if ($null -ne $thermoC) { [Reflection.Assembly]::LoadFrom($thermoC.FullName) | Out-Null }
[Reflection.Assembly]::LoadFrom((Join-Path $engineDirectory "DWSIM.Thermodynamics.dll")) | Out-Null
[Reflection.Assembly]::LoadFrom((Join-Path $engineDirectory "DWSIM.Automation.dll")) | Out-Null
$automation = New-Object DWSIM.Automation.Automation3
$flowsheet = $automation.LoadFlowsheet2($resolvedCase)
$package = [DwsimUnifacExtractor]::PropertyPackage($flowsheet)
[DwsimUnifacExtractor]::SetCurrentMaterialStream($flowsheet, $package)

$basis = @([DwsimUnifacExtractor]::CompoundBasis($package) | ForEach-Object {
    [ordered]@{
        compound_id = $_.Name; q = $_.Q; r = $_.R
        group_surface_fractions = @($_.SurfaceFractions | ForEach-Object {
            [ordered]@{ secondary_id = $_.SecondaryId; value = $_.Value }
        })
    }
})
$groups = @([DwsimUnifacExtractor]::GroupRows($package) | ForEach-Object {
    [ordered]@{
        primary_id = $_.PrimaryId; secondary_id = $_.SecondaryId
        primary_name = $_.PrimaryName; group_name = $_.GroupName; r = $_.R; q = $_.Q
    }
})
$interactions = @([DwsimUnifacExtractor]::InteractionRows($package) | ForEach-Object {
    [ordered]@{ first_primary_id = $_.FirstPrimaryId; second_primary_id = $_.SecondPrimaryId; a_kelvin = $_.A }
})
if ($basis.Count -ne 2 -or $groups.Count -eq 0 -or $interactions.Count -eq 0) {
    throw "UNIFAC extraction returned an invalid record count"
}
$document = [ordered]@{
    schema_version = "dwsim-unifac-data-1"; model = "UNIFAC"
    source = [ordered]@{
        product = "DWSIM"; revision = $DwsimRevision
        case_sha256 = (Get-FileHash -LiteralPath $resolvedCase -Algorithm SHA256).Hash.ToLowerInvariant()
        groups_resource = "DWSIM.Thermodynamics.unifac.txt"
        groups_sha256 = [DwsimUnifacExtractor]::ResourceSha256($package, "DWSIM.Thermodynamics.unifac.txt")
        interactions_resource = "DWSIM.Thermodynamics.unifac_ip.txt"
        interactions_sha256 = [DwsimUnifacExtractor]::ResourceSha256($package, "DWSIM.Thermodynamics.unifac_ip.txt")
    }
    compound_basis = $basis; groups = $groups; interactions = $interactions
}
$resolvedOutput = [IO.Path]::GetFullPath((Join-Path (Get-Location) $OutputPath))
$directory = Split-Path -Parent $resolvedOutput
if (-not (Test-Path -LiteralPath $directory -PathType Container)) {
    New-Item -ItemType Directory -Path $directory -Force | Out-Null
}
$document | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $resolvedOutput -Encoding utf8
[ordered]@{
    output = $resolvedOutput; groups = $groups.Count; directed_interactions = $interactions.Count
    compound_basis_records = $basis.Count
    sha256 = (Get-FileHash -LiteralPath $resolvedOutput -Algorithm SHA256).Hash.ToLowerInvariant()
} | ConvertTo-Json
