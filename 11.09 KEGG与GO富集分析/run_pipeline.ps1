$ErrorActionPreference = 'Stop'
$dataRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$rExecutable = (Get-Command Rscript.exe -ErrorAction Stop).Source
Push-Location -LiteralPath $dataRoot
try {
    python 'enrichment_analysis/00_scripts/download_databases.py'
    if ($LASTEXITCODE -ne 0) { throw 'Database preparation failed' }
    & $rExecutable 'enrichment_analysis/00_scripts/run_enrichment.R'
    if ($LASTEXITCODE -ne 0) { throw 'Enrichment analysis failed' }
    & $rExecutable 'enrichment_analysis/00_scripts/plot_enrichment.R'
    if ($LASTEXITCODE -ne 0) { throw 'Enrichment plotting failed' }
    python 'enrichment_analysis/00_scripts/make_report.py'
    if ($LASTEXITCODE -ne 0) { throw 'Report or validation failed' }
} finally {
    Pop-Location
}
