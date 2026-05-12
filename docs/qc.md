# Quality Control

The collaborator-facing pipeline can generate a QC report automatically:

```matlab
cfg.generateQc = true;
cfg.qcDir = "";          % default: <outputDir>/<prefix>_qc
cfg.motionFile = "";     % optional SPM/FSL/AFNI motion text file
outputs = fmriPhase.pipeline.generateRegressors(cfg);
```

The QC folder contains:

- `qc_summary.txt`: input files, selected SVD component, TR, runtime, and
  optional motion R-squared.
- `regressors.png`: all exported phase regressors.
- `svd_selection.png`: SVD component selection scores.
- `respiratory_phase.png`: selected SVD timecourse and Hilbert respiratory
  phase.
- `breathing_off_resonance_field.png`: low, high, and high-minus-low field map.
- `regressor_correlation.png`: correlation among phase-derived regressors.

If `cfg.motionFile` is provided, the report also includes:

- `motion_correlation.png`: correlation between phase regressors and motion
  model columns.
- `selected_svd_vs_motion.png`: selected SVD timecourse overlaid with an
  approximate framewise-displacement trace.

Motion QC is not meant to reject data automatically. It helps distinguish
smooth breathing-related off-resonance structure from regressors that are mostly
explained by head motion.
