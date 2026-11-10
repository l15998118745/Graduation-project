$ErrorActionPreference = 'Stop'
$dataRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$rExecutable = (Get-Command Rscript.exe -ErrorAction Stop).Source
Push-Location -LiteralPath $dataRoot
try {
    python 'deseq2_baseline_only_analysis/00_scripts/annotate_genes.py'
    if ($LASTEXITCODE -ne 0) { throw 'Gene annotation failed' }
    & $rExecutable 'deseq2_baseline_only_analysis/00_scripts/plot_annotated.R'
    if ($LASTEXITCODE -ne 0) { throw 'Annotated plotting failed' }
    python 'deseq2_baseline_only_analysis/00_scripts/finish_annotation.py'
    if ($LASTEXITCODE -ne 0) { throw 'Final validation failed' }
} finally {
    Pop-Location
}
