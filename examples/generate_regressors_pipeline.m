%% Generate fMRI phase nuisance regressors for one run.
%
% Edit the paths below, then run this script from MATLAB. Outputs include
% generic TSV/MAT files plus SPM, AFNI, and CONN import-friendly files.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'matlab'));
addpath(fullfile(repoRoot, 'matlab', 'vendor', 'tapas_physio'));

cfg = fmriPhase.config();

% Required inputs.
cfg.phaseFile = "/path/to/sub-01_ses-01_task-rest_part-phase_bold.nii";
cfg.magnitudeFile = "/path/to/sub-01_ses-01_task-rest_part-mag_bold.nii";
cfg.qsmSetup = "/path/to/QSM/addpathqsm.m";
cfg.outputDir = "/path/to/phase_regressor_outputs";

% Recommended run identifier. This becomes the output filename prefix.
cfg.regressorPrefix = "sub-01_ses-01_task-rest";

% Acquisition and processing assumptions. Change these to match your data.
cfg.echoTime = 0.030;
cfg.fieldStrength = 3.0;
cfg.gyromagneticRatio = 267.513;
cfg.phaseInputUnits = "auto";
cfg.phaseScale = pi / 2048;
cfg.phaseOffset = -pi;
cfg.rotate90 = false;
cfg.inputImage = "harmonicField";
cfg.interpolationFactor = 1;
cfg.nSvdComponents = 5;
cfg.harmonicOrder = 3;

% Export all supported formats. You can also request a subset, e.g. ["spm"].
cfg.exportFormats = ["matlab", "bids", "spm", "afni", "conn"];

outputs = fmriPhase.pipeline.generateRegressors(cfg);
disp(outputs);
