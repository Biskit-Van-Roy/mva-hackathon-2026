$ErrorActionPreference = "Stop"

$ProjectRoot = Split-Path -Parent $PSScriptRoot
$PythonExe = Join-Path $ProjectRoot ".venv\Scripts\python.exe"
$VcfPath = Join-Path $ProjectRoot "data\raw\WGS_EX2312012_HGWCNDSX7.vcf.gz"
$PhenotypePath = Join-Path $ProjectRoot "data\raw\Challenge_Clinical_Phenotype_1.docx"
$OutputPath = Join-Path $ProjectRoot "artifacts\intake_report.json"

if (-not (Test-Path $PythonExe)) {
    throw "The environment is missing. Run ./scripts/setup.ps1 first."
}

& $PythonExe -m mva_hackathon.intake `
    --vcf $VcfPath `
    --phenotype $PhenotypePath `
    --output $OutputPath

if ($LASTEXITCODE -ne 0) {
    throw "Input inspection failed."
}
