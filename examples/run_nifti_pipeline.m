%% Example MATLAB v2 fMRI phase pipeline.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'matlab'));
addpath(fullfile(repoRoot, 'matlab', 'vendor', 'tapas_physio'));

cfg = fmriPhase.config();
cfg.phaseFile = "/path/to/sub-01_task-rest_part-phase_bold.nii";
cfg.magnitudeFile = "/path/to/sub-01_task-rest_part-mag_bold.nii";
cfg.outputDir = fullfile(repoRoot, "outputs", "sub-01");
cfg.qsmSetup = "/path/to/QSM/addpathqsm.m";
cfg.echoTime = 0.030;
cfg.fieldStrength = 3.0;
cfg.rotate90 = false;
cfg.interpolationFactor = 1;
cfg.doPlot = false;

result = fmriPhase.pipeline.run(cfg);
disp(size(result.regressors));
