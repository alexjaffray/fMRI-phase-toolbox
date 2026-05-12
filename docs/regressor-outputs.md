# Regressor Outputs

Use `fmriPhase.pipeline.generateRegressors(cfg)` to generate fMRI
phase-derived nuisance regressors for a single run.

The primary matrix contains:

- `phase_harmonic_01` ... `phase_harmonic_16`: solid-harmonic coefficients
  fitted to the respiratory-dominant phase field.
- `phase_svd_01` ... `phase_svd_05`: leading SVD component timecourses used
  during respiratory component selection.

All columns are normalized with MATLAB `normalize(..., 1)`, matching the legacy
scripts.

## Files

For a prefix such as `sub-01_ses-01_task-rest`, the exporter writes:

- `*_phase_regressors.mat`: MATLAB archive with `regressors`, `R`, names, TR,
  and config.
- `*_desc-fMRIPhase_timeseries.tsv`: tab-separated matrix with column headers.
- `*_desc-fMRIPhase_timeseries.json`: metadata sidecar for the TSV.
- `*_spm_multiple_regressors.mat`: SPM-style MAT file containing `R`.
- `*_spm_multiple_regressors.txt`: numeric tab-delimited matrix for SPM.
- `*_afni.1D`: AFNI 1D matrix with column labels in comments.
- `*_conn_covariates.tsv`: multi-column table for CONN import workflows.
- `*_conn_covariates.mat`: MATLAB copy with `data` and `namesCell`.
- `*_conn_covariates/*.txt`: one text file per covariate, useful when importing
  covariates one at a time.

## Minimal Example

```matlab
addpath('/path/to/fMRI-phase-toolbox-v2/matlab');

cfg = fmriPhase.config();
cfg.phaseFile = "/data/sub-01_task-rest_part-phase_bold.nii";
cfg.magnitudeFile = "/data/sub-01_task-rest_part-mag_bold.nii";
cfg.qsmSetup = "/opt/QSM/addpathqsm.m";
cfg.outputDir = "/data/derivatives/fmri_phase/sub-01";
cfg.regressorPrefix = "sub-01_task-rest";

outputs = fmriPhase.pipeline.generateRegressors(cfg);
```
