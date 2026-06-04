function cfg = validateConfig(cfg)
%VALIDATECONFIG Validate and normalize pipeline configuration.

required = ["phaseFile", "magnitudeFile"];
for name = required
    if strlength(string(cfg.(name))) == 0
        error('fmriPhase:config:MissingRequired', ...
            'Configuration field "%s" is required.', name);
    end
end

if cfg.dataFormat ~= "nifti"
    error('fmriPhase:config:UnsupportedDataFormat', ...
        'Only dataFormat="nifti" is implemented in v2 so far.');
end

cfg.phaseInputUnits = string(cfg.phaseInputUnits);
validPhaseInputUnits = ["auto", "radians", "scaled"];
if ~any(cfg.phaseInputUnits == validPhaseInputUnits)
    error('fmriPhase:config:PhaseInputUnits', ...
        'phaseInputUnits must be "auto", "radians", or "scaled".');
end

if ~(isempty(cfg.phaseScale) || isFiniteNumericScalar(cfg.phaseScale))
    error('fmriPhase:config:PhaseScale', ...
        'phaseScale must be empty or a finite numeric scalar.');
end

if cfg.phaseInputUnits == "scaled" && isempty(cfg.phaseScale)
    error('fmriPhase:config:PhaseScaleRequired', ...
        'phaseScale must be specified when phaseInputUnits="scaled".');
end

if ~isFiniteNumericScalar(cfg.phaseOffset)
    error('fmriPhase:config:PhaseOffset', ...
        'phaseOffset must be a finite numeric scalar.');
end

if ~isFiniteNumericScalar(cfg.phaseRangeTolerance) || cfg.phaseRangeTolerance < 0
    error('fmriPhase:config:PhaseRangeTolerance', ...
        'phaseRangeTolerance must be a nonnegative finite numeric scalar.');
end

if cfg.interpolationFactor < 1 || fix(cfg.interpolationFactor) ~= cfg.interpolationFactor
    error('fmriPhase:config:InterpolationFactor', ...
        'interpolationFactor must be a positive integer.');
end

if cfg.nSvdComponents < 1 || fix(cfg.nSvdComponents) ~= cfg.nSvdComponents
    error('fmriPhase:config:SvdComponents', ...
        'nSvdComponents must be a positive integer.');
end

if strlength(string(cfg.outputDir)) == 0
    cfg.outputDir = pwd;
end

if ~isfolder(cfg.outputDir)
    mkdir(cfg.outputDir);
end

end

function tf = isFiniteNumericScalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value);
end
