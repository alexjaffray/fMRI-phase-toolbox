# fMRI Phase Toolbox 

I've been a bit lazy in updating this package, but this is a v0.2 release of the fMRI phase toolbox for processing spatiotemporal image phase fluctuations during fMRI scans. It contains an implementation of the method described in:

Detection of respiration-induced field modulations in fMRI: A concurrent and navigator-free approach https://doi.org/10.1162/imag_a_00091 

This v0.2 release cleans up a lot of older code and research scripts, and provides a configuration interface for setting up the method. I have been wanting to do this for a while, and also had been hearing about how good coding agents have gotten recently, so I decided to use this as a chance to both clean up the repo and try to see how working together with the machine goes. This is a work in progress, all feedback is welcome.

## List of Changes against Initial Release

- MATLAB package namespace: `+fmriPhase`
- Non-interactive configuration struct instead of hard-coded local paths
- Export of package-friendly nuisance regressors for SPM, AFNI, and CONN
- Automatic QC figures, progress messages, and optional motion comparison
- Legacy MATLAB source copied under `legacy/` for comparison

## Core Features

- QSM-backed phase unwrapping and background field removal
- SVD-based respiratory component extraction
- Spherical-harmonic regressor generation

## Dependencies

Required for the current MATLAB pipeline:

- MATLAB with Image Processing Toolbox and Signal Processing Toolbox
- Christian Kames' QSM toolbox (`QSM.m`) for `generateMask`,
  `unwrapLaplacian`, and `resharp`
- NIfTI support through MATLAB's `niftiread` and `niftiinfo`
- FSL Brain Extraction Tool accessible via `bet` command on CLI

## MATLAB Quick Start

```matlab
addpath('/path/to/this/toolbox/matlab');
addpath('/path/to/this/toolbox/matlab/vendor/tapas_physio');

cfg = fmriPhase.config();
cfg.phaseFile = "/path/to/sub-01_part-phase_bold.nii";
cfg.magnitudeFile = "/path/to/sub-01_part-mag_bold.nii";
cfg.outputDir = "/path/to/output";
cfg.qsmSetup = "/path/to/QSM/addpathqsm.m";
cfg.generateQc = true;

% Phase input handling. "auto" treats images with range near 2*pi as radians;
% otherwise it applies phaseRadians = phaseImage * phaseScale + phaseOffset.
% Defaults preserve the legacy integer NIfTI scaling: pi / 2048, then -pi.
cfg.phaseInputUnits = "auto"; % "auto", "radians", or "scaled"
cfg.phaseScale = pi / 2048;
cfg.phaseOffset = -pi;

result = fmriPhase.pipeline.run(cfg);
```

Export of regressors is available so that they may be used or shared with collaborators. My research primarily focuses on MRI methods, so I'm quite open to suggestions on how to make this more useable for the neuroscience community.

```matlab
cfg.regressorPrefix = "sub-01_task-rest";
outputs = fmriPhase.pipeline.generateRegressors(cfg);
```

See `examples/run_nifti_pipeline.m` and
`examples/generate_regressors_collaborator.m` for fuller examples.


## Layout

- `matlab/+fmriPhase/`: MATLAB package
- `matlab/vendor/tapas_physio/`: copied TAPAS/PhysIO-derived helpers retained
  with their upstream headers
- `legacy/`: reference copy of the original MATLAB source and pipeline scripts
- `docs/`: refactor notes and roadmap
- `tests/`: MATLAB unit tests for package-level functions

## Note

TAPAS/PhysIO-derived files are included for consistency and research workflow, and can be found in `matlab/vendor/tapas_physio/`, for completeness. However, users are encouraged to use the PhysIO toolbox tools that best fit their research needs, in the spirit of TAPAS.

The original repository maintained by Lars Kasper (@mrikasper) can be found at: https://github.com/ComputationalPsychiatry/PhysIO 
