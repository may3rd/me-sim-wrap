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

public sealed class UniquacSourceRow
{
    public int FirstId { get; set; }
    public int SecondId { get; set; }
    public double A12 { get; set; }
    public double A21 { get; set; }
    public string Comment { get; set; }
}

public sealed class UniquacCompoundBasis
{
    public string Name { get; set; }
    public string ChemSepId { get; set; }
    public double Q { get; set; }
    public double R { get; set; }
}

public sealed class UniquacResolvedPair
{
    public string First { get; set; }
    public string Second { get; set; }
    public double A12 { get; set; }
    public double A21 { get; set; }
    public double B12 { get; set; }
    public double B21 { get; set; }
    public double C12 { get; set; }
    public double C21 { get; set; }
}

public static class DwsimUniquacExtractor
{
    private const BindingFlags PublicInstance = BindingFlags.Public | BindingFlags.Instance;
    private const string ResourceName = "DWSIM.Thermodynamics.uniquac.dat";

    private static object Get(object target, string name)
    {
        PropertyInfo property = target.GetType().GetProperties(PublicInstance)
            .FirstOrDefault(candidate =>
                candidate.Name == name && candidate.GetIndexParameters().Length == 0);
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
            if (String.Equals(
                    Convert.ToString(Get(entry.Value, "ComponentName"), CultureInfo.InvariantCulture),
                    "UNIQUAC", StringComparison.Ordinal)) return entry.Value;
        throw new InvalidOperationException("No exact UNIQUAC property-package match was found.");
    }

    public static void SetCurrentMaterialStream(object flowsheet, object propertyPackage)
    {
        object stream = null;
        foreach (DictionaryEntry entry in (IDictionary)Get(flowsheet, "SimulationObjects"))
            if (entry.Value != null && entry.Value.GetType().FullName ==
                    "DWSIM.Thermodynamics.Streams.MaterialStream") stream = entry.Value;
        if (stream == null) throw new InvalidOperationException("The case has no material stream.");
        propertyPackage.GetType().GetProperties(PublicInstance)
            .First(candidate => candidate.Name == "CurrentMaterialStream" && candidate.CanWrite)
            .SetValue(propertyPackage, stream, null);
    }

    private static double[] Values(object package, string method)
    {
        return ((IEnumerable)Invoke(package, method)).Cast<object>()
            .Select(value => Convert.ToDouble(value, CultureInfo.InvariantCulture)).ToArray();
    }

    private static string[] Strings(object package, string method)
    {
        return ((IEnumerable)Invoke(package, method)).Cast<object>()
            .Select(value => Convert.ToString(value, CultureInfo.InvariantCulture)).ToArray();
    }

    public static UniquacCompoundBasis[] CompoundBasis(object package)
    {
        string[] names = Strings(package, "RET_VNAMES");
        string[] ids = Strings(package, "RET_VIDS");
        double[] q = Values(package, "RET_VQ");
        double[] r = Values(package, "RET_VR");
        return names.Select((name, index) => new UniquacCompoundBasis {
            Name = name, ChemSepId = ids[index], Q = q[index], R = r[index]
        }).ToArray();
    }

    public static UniquacResolvedPair ResolvedPair(object package)
    {
        string[] names = Strings(package, "RET_VNAMES");
        object model = Invoke(package, "GetModel");
        IDictionary outer = (IDictionary)Get(model, "InteractionParameters");
        object pair = ((IDictionary)outer[names[0]])[names[1]];
        return new UniquacResolvedPair {
            First = names[0], Second = names[1],
            A12 = Convert.ToDouble(Get(pair, "A12"), CultureInfo.InvariantCulture),
            A21 = Convert.ToDouble(Get(pair, "A21"), CultureInfo.InvariantCulture),
            B12 = Convert.ToDouble(Get(pair, "B12"), CultureInfo.InvariantCulture),
            B21 = Convert.ToDouble(Get(pair, "B21"), CultureInfo.InvariantCulture),
            C12 = Convert.ToDouble(Get(pair, "C12"), CultureInfo.InvariantCulture),
            C21 = Convert.ToDouble(Get(pair, "C21"), CultureInfo.InvariantCulture)
        };
    }

    public static UniquacSourceRow[] SourceRows(object package)
    {
        var rows = new List<UniquacSourceRow>();
        using (Stream stream = package.GetType().Assembly.GetManifestResourceStream(ResourceName))
        using (var reader = new StreamReader(stream))
            while (!reader.EndOfStream)
            {
                string line = reader.ReadLine();
                if (String.IsNullOrWhiteSpace(line)) continue;
                string[] fields = line.Split(new[] {';'}, 5);
                int firstId;
                int secondId;
                if (!Int32.TryParse(fields[0], NumberStyles.Integer,
                        CultureInfo.InvariantCulture, out firstId) ||
                    !Int32.TryParse(fields[1], NumberStyles.Integer,
                        CultureInfo.InvariantCulture, out secondId)) continue;
                rows.Add(new UniquacSourceRow {
                    FirstId = firstId,
                    SecondId = secondId,
                    A12 = Double.Parse(fields[2], CultureInfo.InvariantCulture),
                    A21 = Double.Parse(fields[3], CultureInfo.InvariantCulture),
                    Comment = fields.Length == 5 ? fields[4] : ""
                });
            }
        return rows.ToArray();
    }

    public static string ResourceSha256(object package)
    {
        using (Stream stream = package.GetType().Assembly.GetManifestResourceStream(ResourceName))
        using (SHA256 sha = SHA256.Create())
            return String.Concat(sha.ComputeHash(stream).Select(value => value.ToString("x2")));
    }
}
"@

$engineDirectory = (Resolve-Path -LiteralPath $EngineBin).Path
$resolvedCase = (Resolve-Path -LiteralPath $CasePath).Path
foreach ($file in @("DWSIM.Interfaces.dll", "DWSIM.Thermodynamics.dll", "DWSIM.Automation.dll")) {
    $path = Join-Path $engineDirectory $file
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Missing DWSIM assembly: $path"
    }
}
[Reflection.Assembly]::LoadFrom((Join-Path $engineDirectory "DWSIM.Interfaces.dll")) | Out-Null
$thermoC = Get-ChildItem -Path $engineDirectory -Filter "ThermoCS.dll" -File -Recurse `
    -ErrorAction SilentlyContinue | Select-Object -First 1
if ($null -ne $thermoC) { [Reflection.Assembly]::LoadFrom($thermoC.FullName) | Out-Null }
[Reflection.Assembly]::LoadFrom((Join-Path $engineDirectory "DWSIM.Thermodynamics.dll")) | Out-Null
[Reflection.Assembly]::LoadFrom((Join-Path $engineDirectory "DWSIM.Automation.dll")) | Out-Null

$automation = New-Object DWSIM.Automation.Automation3
$flowsheet = $automation.LoadFlowsheet2($resolvedCase)
$package = [DwsimUniquacExtractor]::PropertyPackage($flowsheet)
[DwsimUniquacExtractor]::SetCurrentMaterialStream($flowsheet, $package)
$basis = @([DwsimUniquacExtractor]::CompoundBasis($package) | ForEach-Object {
    [ordered]@{ compound_id = $_.Name; chemsep_id = $_.ChemSepId; q = $_.Q; r = $_.R }
})
$pair = [DwsimUniquacExtractor]::ResolvedPair($package)
$rows = @([DwsimUniquacExtractor]::SourceRows($package) | ForEach-Object {
    [ordered]@{
        first_chemsep_id = $_.FirstId; second_chemsep_id = $_.SecondId
        A12 = $_.A12; A21 = $_.A21; comment = $_.Comment
    }
})
if ($basis.Count -ne 2 -or $rows.Count -eq 0) {
    throw "UNIQUAC extraction returned an invalid basis or source-table count"
}

$document = [ordered]@{
    schema_version = "dwsim-uniquac-data-1"
    model = "UNIQUAC"
    source = [ordered]@{
        product = "DWSIM"; revision = $DwsimRevision
        case_sha256 = (Get-FileHash -LiteralPath $resolvedCase -Algorithm SHA256).Hash.ToLowerInvariant()
        resource = "DWSIM.Thermodynamics.uniquac.dat"
        resource_sha256 = [DwsimUniquacExtractor]::ResourceSha256($package)
        interaction_basis = "cal/mol"; gas_constant_cal_per_mol_k = 1.98721
    }
    compound_basis = $basis
    resolved_interaction = [ordered]@{
        first = $pair.First; second = $pair.Second
        A12 = $pair.A12; A21 = $pair.A21
        B12 = $pair.B12; B21 = $pair.B21; C12 = $pair.C12; C21 = $pair.C21
    }
    source_interactions = $rows
}
$resolvedOutput = [IO.Path]::GetFullPath((Join-Path (Get-Location) $OutputPath))
$directory = Split-Path -Parent $resolvedOutput
if (-not (Test-Path -LiteralPath $directory -PathType Container)) {
    New-Item -ItemType Directory -Path $directory -Force | Out-Null
}
$document | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $resolvedOutput -Encoding utf8
[ordered]@{
    output = $resolvedOutput; source_pairs = $rows.Count; compound_basis_records = $basis.Count
    sha256 = (Get-FileHash -LiteralPath $resolvedOutput -Algorithm SHA256).Hash.ToLowerInvariant()
} | ConvertTo-Json
