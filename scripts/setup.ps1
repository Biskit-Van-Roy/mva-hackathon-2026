$ErrorActionPreference = "Stop"

$ProjectRoot = Split-Path -Parent $PSScriptRoot
Set-Location $ProjectRoot

if (Get-Command py -ErrorAction SilentlyContinue) {
    & py -3 -m venv .venv
}
elseif (Get-Command python -ErrorAction SilentlyContinue) {
    & python -m venv .venv
}
else {
    throw "Python 3.11 or 3.12 was not found in PATH."
}

$PythonExe = Join-Path $ProjectRoot ".venv\Scripts\python.exe"
& $PythonExe -m pip install --upgrade pip
& $PythonExe -m pip install -e ".[download]"

Write-Host "Environment ready. Next: ./.venv/Scripts/hf.exe auth login"
