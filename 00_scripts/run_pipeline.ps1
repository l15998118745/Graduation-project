$ErrorActionPreference = 'Stop'
$analysisRoot = Split-Path -Parent $PSScriptRoot
$dataRoot = Split-Path -Parent $analysisRoot
$rExecutable = (Get-Command Rscript.exe -ErrorAction Stop).Source
Push-Location -LiteralPath $dataRoot
try {
    & $rExecutable 'deseq2_baseline_only_analysis/00_scripts/run_deseq2.R'
    if ($LASTEXITCODE -ne 0) { throw 'DESeq2 pipeline failed; inspect the R error.' }
    & $rExecutable 'deseq2_baseline_only_analysis/00_scripts/summarize_results.R'
    if ($LASTEXITCODE -ne 0) { throw 'Result summarization failed; inspect the R error.' }
    python 'deseq2_baseline_only_analysis/00_scripts/make_report.py'
    if ($LASTEXITCODE -ne 0) { throw 'Markdown report failed; check Python and pandas.' }
    python 'deseq2_baseline_only_analysis/00_scripts/organize_outputs.py' --source 'deseq2_baseline_only_analysis/.work'
    if ($LASTEXITCODE -ne 0) { throw 'Output organization failed' }
    & "$PSScriptRoot/refresh_annotations.ps1"
} finally {
    Pop-Location
}
