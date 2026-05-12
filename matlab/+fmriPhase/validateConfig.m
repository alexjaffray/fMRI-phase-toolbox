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
