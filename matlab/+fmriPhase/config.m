function cfg = config(varargin)
%CONFIG Create a default configuration for the fMRI phase pipeline.
%
% cfg = fmriPhase.config() returns defaults.
% cfg = fmriPhase.config('phaseFile', "...") overrides named fields.

cfg = struct();
cfg.phaseFile = "";
cfg.magnitudeFile = "";
cfg.physLogFile = "";
cfg.outputDir = "";
cfg.qsmSetup = "";
cfg.mreconSetup = "";
cfg.dataFormat = "nifti";
cfg.rotate90 = false;
cfg.echoTime = 0.030;
cfg.fieldStrength = 3.0;
cfg.gyromagneticRatio = 267.513;
cfg.b0Direction = [0 0 1];
cfg.inputImage = "harmonicField";
cfg.interpolationFactor = 1;
cfg.nSvdComponents = 5;
cfg.respComponentMethod = "maxDerivativeL1";
cfg.doPlot = false;
cfg.maskOptions = '-m -n -f 0.5';
cfg.resharpRadii = [];
cfg.resharpRegularization = 0.05;
cfg.harmonicOrder = 3;
cfg.saveRegressors = true;
cfg.regressorFilename = "fmri_phase_regressors.mat";
cfg.regressorPrefix = "";
cfg.exportFormats = ["matlab", "bids", "spm", "afni", "conn"];
cfg.verbose = true;
cfg.generateQc = true;
cfg.qcDir = "";
cfg.motionFile = "";
cfg.motionModel = "motion+derivatives";
cfg.physioFile = "";

if mod(numel(varargin), 2) ~= 0
    error('fmriPhase:config:NameValuePairs', ...
        'Overrides must be provided as name-value pairs.');
end

for idx = 1:2:numel(varargin)
    name = varargin{idx};
    value = varargin{idx + 1};
    if ~isfield(cfg, name)
        error('fmriPhase:config:UnknownField', ...
            'Unknown configuration field "%s".', name);
    end
    cfg.(name) = value;
end

end
