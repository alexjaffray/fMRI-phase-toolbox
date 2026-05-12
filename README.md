# fMRI Phase Toolbox 

Full refactor of the fMRI phase toolbox for processing spatiotemporal
image phase fluctuations during fMRI scans.

The original project at `/srv/data/ajaffray/fMRI-phase-toolbox` is preserved
intact. This v2 folder separates reusable toolbox code from legacy research
scripts, data, and paper outputs.

## Current Scope

- MATLAB package namespace: `+fmriPhase`
- Non-interactive configuration struct instead of hard-coded local paths
- NIfTI phase/magnitude loading
- QSM-backed phase unwrapping and background field removal
- SVD-based respiratory component extraction
- Spherical-harmonic regressor generation
- Export of package-friendly nuisance regressors for SPM, AFNI, and CONN
- Automatic QC figures, progress messages, and optional motion comparison
- Legacy MATLAB source copied under `legacy/` for comparison

## MATLAB Quick Start

```matlab
addpath('/srv/data/ajaffray/fMRI-phase-toolbox-v2/matlab');
addpath('/srv/data/ajaffray/fMRI-phase-toolbox-v2/matlab/vendor/tapas_physio');

cfg = fmriPhase.config();
cfg.phaseFile = "/path/to/sub-01_part-phase_bold.nii";
cfg.magnitudeFile = "/path/to/sub-01_part-mag_bold.nii";
cfg.outputDir = "/path/to/output";
cfg.qsmSetup = "/path/to/QSM/addpathqsm.m";
cfg.generateQc = true;

result = fmriPhase.pipeline.run(cfg);
```

For collaborator-facing regressor export, use:

```matlab
cfg.regressorPrefix = "sub-01_task-rest";
outputs = fmriPhase.pipeline.generateRegressors(cfg);
```

See `examples/run_nifti_pipeline.m` and
`examples/generate_regressors_collaborator.m` for fuller examples.

## Dependencies

Required for the current MATLAB pipeline:

- MATLAB with Image Processing Toolbox and Signal Processing Toolbox
- Christian Kames' QSM toolbox (`QSM.m`) for `generateMask`,
  `unwrapLaplacian`, and `resharp`
- NIfTI support through MATLAB's `niftiread` and `niftiinfo`

Optional:

- Julia/Python ports are not started yet; the first target is a stable MATLAB
  API with testable numerical boundaries

## Layout

- `matlab/+fmriPhase/`: MATLAB package
- `matlab/vendor/tapas_physio/`: copied TAPAS/PhysIO-derived helpers retained
  with their upstream headers
- `legacy/`: reference copy of the original MATLAB source and pipeline scripts
- `docs/`: refactor notes and roadmap
- `tests/`: MATLAB unit tests for package-level functions

## Licensing Note

The original project is MIT licensed. Some TAPAS/PhysIO-derived files included
in the original source state GPL terms in their headers; v2 keeps those files in
`matlab/vendor/tapas_physio/` so downstream licensing decisions are explicit.
