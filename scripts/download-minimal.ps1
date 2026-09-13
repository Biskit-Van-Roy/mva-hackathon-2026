param(
    [string]$Token = $env:HF_TOKEN
)

$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$HfExe = Join-Path $ProjectRoot ".venv\Scripts\hf.exe"
$Destination = Join-Path $ProjectRoot "data\raw"
$Repository = "SageBio/mva-hackathon-2026-data"

if (-not (Test-Path $HfExe)) {
    throw "The environment is missing. Run ./scripts/setup.ps1 first."
}

$Files = @(
    "Challenge_Clinical_Phenotype_1.docx",
    "WGS_EX2312012_HGWCNDSX7.vcf.gz",
    "WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi"
)

foreach ($File in $Files) {
    $Arguments = @(
        "download",
        $Repository,
        $File,
        "--repo-type", "dataset",
        "--local-dir", $Destination
    )
    if ($Token) {
        $Arguments += @("--token", $Token)
    }
    & $HfExe @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Download failed for $File. Confirm that dataset access was approved."
    }
}

Write-Host "Minimal data downloaded to $Destination"
