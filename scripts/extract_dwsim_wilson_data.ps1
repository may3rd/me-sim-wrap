[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$EngineBin,

    [Parameter(Mandatory = $true)]
    [string]$CasePath,

    [Parameter(Mandatory = $true)]
    [string]$OutputPath,

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

public sealed class WilsonInteractionRow
{
    public string FirstCas { get; set; }
    public string SecondCas { get; set; }
    public double A12 { get; set; }
    public double A21 { get; set; }
}

public sealed class WilsonCompoundBasis
{
    public string Name { get; set; }
    public string Cas { get; set; }
    public double MolarVolume { get; set; }
}

public static class DwsimWilsonExtractor
{
    private const BindingFlags PublicInstance = BindingFlags.Public | BindingFlags.Instance;

    private static object Get(object target, string name)
    {
        PropertyInfo property = target.GetType().GetProperties(PublicInstance)
            .FirstOrDefault(candidate =>
                candidate.Name == name && candidate.GetIndexParameters().Length == 0);
        if (property != null) return property.GetValue(target, null);
        FieldInfo field = target.GetType().GetField(name, PublicInstance);
        return field == null ? null : field.GetValue(target);
    }

    private static object Invoke(object target, string name, object[] arguments, int count)
    {
        MethodInfo method = target.GetType().GetMethods(PublicInstance)
            .First(candidate => candidate.Name == name && candidate.GetParameters().Length == count);
        return method.Invoke(target, arguments);
    }

    public static object PropertyPackage(object flowsheet, string componentName)
    {
        IDictionary packages = (IDictionary)Get(flowsheet, "PropertyPackages");
        foreach (DictionaryEntry entry in packages)
        {
            object value = entry.Value;
            if (String.Equals(
                    Convert.ToString(Get(value, "ComponentName"), CultureInfo.InvariantCulture),
                    componentName,
                    StringComparison.Ordinal)) return value;
        }
        throw new InvalidOperationException("No exact Wilson property-package match was found.");
    }

    public static void SetCurrentMaterialStream(object flowsheet, object propertyPackage)
    {
        IDictionary objects = (IDictionary)Get(flowsheet, "SimulationObjects");
        object stream = null;
        foreach (DictionaryEntry entry in objects)
        {
            if (entry.Value != null && entry.Value.GetType().FullName ==
                    "DWSIM.Thermodynamics.Streams.MaterialStream")
            {
                stream = entry.Value;
                break;
            }
        }
        if (stream == null) throw new InvalidOperationException("The case has no material stream.");
        PropertyInfo property = propertyPackage.GetType().GetProperties(PublicInstance)
            .First(candidate => candidate.Name == "CurrentMaterialStream" && candidate.CanWrite);
        property.SetValue(propertyPackage, stream, null);
    }

    public static WilsonInteractionRow[] Interactions(object propertyPackage)
    {
        object model = Get(propertyPackage, "WilsonM");
        IDictionary outer = (IDictionary)Get(model, "BIPs");
        var rows = new List<WilsonInteractionRow>();
        foreach (DictionaryEntry first in outer)
        {
            IDictionary inner = (IDictionary)first.Value;
            foreach (DictionaryEntry second in inner)
            {
                double[] values = (double[])second.Value;
                rows.Add(new WilsonInteractionRow {
                    FirstCas = (string)first.Key,
                    SecondCas = (string)second.Key,
                    A12 = values[0],
                    A21 = values[1]
                });
            }
        }
        return rows.OrderBy(row => row.FirstCas, StringComparer.Ordinal)
            .ThenBy(row => row.SecondCas, StringComparer.Ordinal).ToArray();
    }

    public static WilsonCompoundBasis[] CompoundBasis(object propertyPackage)
    {
        string[] names = ((IEnumerable)Invoke(propertyPackage, "RET_VNAMES", new object[0], 0))
            .Cast<object>().Select(value => Convert.ToString(value, CultureInfo.InvariantCulture)).ToArray();
        string[] cases = ((IEnumerable)Invoke(propertyPackage, "RET_VCAS", new object[0], 0))
            .Cast<object>().Select(value => Convert.ToString(value, CultureInfo.InvariantCulture)).ToArray();
        object[] arguments = ((IEnumerable)Invoke(propertyPackage, "GetArguments", new object[0], 0))
            .Cast<object>().ToArray();
        double[] volumes = ((IEnumerable)arguments[1]).Cast<object>()
            .Select(value => Convert.ToDouble(value, CultureInfo.InvariantCulture)).ToArray();
        return names.Select((name, index) => new WilsonCompoundBasis {
            Name = name,
            Cas = cases[index],
            MolarVolume = volumes[index]
        }).ToArray();
    }

    public static string ResourceSha256(object propertyPackage)
    {
        Assembly assembly = propertyPackage.GetType().Assembly;
        using (Stream stream = assembly.GetManifestResourceStream(
            "DWSIM.Thermodynamics.wilson_bips.csv"))
        using (SHA256 sha = SHA256.Create())
            return String.Concat(sha.ComputeHash(stream).Select(value => value.ToString("x2")));
    }
}
"@

$engineDirectory = (Resolve-Path -LiteralPath $EngineBin).Path
$resolvedCase = (Resolve-Path -LiteralPath $CasePath).Path
$interfacesPath = Join-Path $engineDirectory "DWSIM.Interfaces.dll"
$thermodynamicsPath = Join-Path $engineDirectory "DWSIM.Thermodynamics.dll"
$automationPath = Join-Path $engineDirectory "DWSIM.Automation.dll"
foreach ($assemblyPath in @($interfacesPath, $thermodynamicsPath, $automationPath)) {
    if (-not (Test-Path -LiteralPath $assemblyPath -PathType Leaf)) {
        throw "Missing DWSIM assembly: $assemblyPath"
    }
}

[Reflection.Assembly]::LoadFrom($interfacesPath) | Out-Null
$thermoCAssemblyPath = Get-ChildItem `
    -Path $engineDirectory `
    -Filter "ThermoCS.dll" `
    -File `
    -Recurse `
    -ErrorAction SilentlyContinue |
    Select-Object -First 1
if ($null -ne $thermoCAssemblyPath) {
    [Reflection.Assembly]::LoadFrom($thermoCAssemblyPath.FullName) | Out-Null
}
[Reflection.Assembly]::LoadFrom($thermodynamicsPath) | Out-Null
[Reflection.Assembly]::LoadFrom($automationPath) | Out-Null

$automation = New-Object DWSIM.Automation.Automation3
$flowsheet = $automation.LoadFlowsheet2($resolvedCase)
$propertyPackage = [DwsimWilsonExtractor]::PropertyPackage($flowsheet, "Wilson")
[DwsimWilsonExtractor]::SetCurrentMaterialStream($flowsheet, $propertyPackage)

$interactions = @(
    [DwsimWilsonExtractor]::Interactions($propertyPackage) | ForEach-Object {
        [ordered]@{
            first_cas = $_.FirstCas
            second_cas = $_.SecondCas
            A12 = [ordered]@{ value = $_.A12; unit = "cal/mol" }
            A21 = [ordered]@{ value = $_.A21; unit = "cal/mol" }
        }
    }
)
$compoundBasis = @(
    [DwsimWilsonExtractor]::CompoundBasis($propertyPackage) | ForEach-Object {
        [ordered]@{
            compound_id = $_.Name
            cas = $_.Cas
            molar_volume_298_15_k = [ordered]@{
                value = $_.MolarVolume
                unit = "m3/kmol"
            }
        }
    }
)
if ($interactions.Count -eq 0 -or $compoundBasis.Count -ne 2) {
    throw "Wilson extraction returned an invalid interaction or compound-basis count"
}

$document = [ordered]@{
    schema_version = "dwsim-wilson-data-1"
    model = "Wilson"
    source = [ordered]@{
        product = "DWSIM"
        revision = $DwsimRevision
        case_sha256 = (Get-FileHash -LiteralPath $resolvedCase -Algorithm SHA256).Hash.ToLowerInvariant()
        resource = "DWSIM.Thermodynamics.wilson_bips.csv"
        resource_sha256 = [DwsimWilsonExtractor]::ResourceSha256($propertyPackage)
        interaction_basis = "cal/mol"
        molar_volume_temperature_k = 298.15
    }
    compound_basis = $compoundBasis
    interactions = $interactions
}

$resolvedOutput = [IO.Path]::GetFullPath((Join-Path (Get-Location) $OutputPath))
$outputDirectory = Split-Path -Parent $resolvedOutput
if (-not (Test-Path -LiteralPath $outputDirectory -PathType Container)) {
    New-Item -ItemType Directory -Path $outputDirectory -Force | Out-Null
}
$document | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $resolvedOutput -Encoding utf8

[ordered]@{
    output = $resolvedOutput
    interaction_pairs = $interactions.Count
    compound_basis_records = $compoundBasis.Count
    sha256 = (Get-FileHash -LiteralPath $resolvedOutput -Algorithm SHA256).Hash.ToLowerInvariant()
} | ConvertTo-Json
