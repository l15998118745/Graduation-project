$ErrorActionPreference = 'Stop'
$scriptDir = $PSScriptRoot
$analysisRoot = Split-Path -Parent $scriptDir
$dataRoot = Split-Path -Parent $analysisRoot
$rExecutable = (Get-Command Rscript.exe -ErrorAction Stop).Source
Push-Location -LiteralPath $dataRoot
try {
    python 'deseq2_baseline_only_analysis/00_scripts/fetch_kegg.py'
    if ($LASTEXITCODE -ne 0) { throw 'KEGG cache download failed' }
    & $rExecutable 'deseq2_baseline_only_analysis/00_scripts/run_enrichment.R'
    if ($LASTEXITCODE -ne 0) { throw 'GO/KEGG ORA failed; inspect enrichment_run.log' }
    python 'deseq2_baseline_only_analysis/00_scripts/summarize_enrichment.py'
    if ($LASTEXITCODE -ne 0) { throw 'Enrichment summary validation failed' }
} finally { Pop-Location }
