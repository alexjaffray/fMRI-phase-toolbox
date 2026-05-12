function result = run(cfg)
%RUN Run the MATLAB fMRI phase regressor pipeline.

pipelineTimer = tic;
cfg = fmriPhase.validateConfig(cfg);
fmriPhase.util.log(cfg, 'Starting fMRI phase regressor pipeline.');

if strlength(string(cfg.qsmSetup)) > 0
    stepTimer = tic;
    fmriPhase.util.log(cfg, 'Adding QSM toolbox paths.');
    run(cfg.qsmSetup);
    fmriPhase.util.logElapsed(cfg, 'QSM setup', stepTimer);
end

stepTimer = tic;
fmriPhase.util.log(cfg, 'Loading phase and magnitude NIfTI files.');
data = fmriPhase.io.loadNiftiPair(cfg.phaseFile, cfg.magnitudeFile, cfg.rotate90);
fmriPhase.util.logElapsed(cfg, 'NIfTI loading', stepTimer);

stepTimer = tic;
fmriPhase.util.log(cfg, 'Generating brain mask from final magnitude volume.');
mask0 = generateMask(data.magnitude(:, :, :, end), data.voxelSize, cfg.maskOptions);
fmriPhase.util.logElapsed(cfg, 'Mask generation', stepTimer);

fieldParams = struct( ...
    'echoTime', cfg.echoTime, ...
    'fieldStrength', cfg.fieldStrength, ...
    'gyromagneticRatio', cfg.gyromagneticRatio, ...
    'voxelSize', data.voxelSize, ...
    'resharpRadii', cfg.resharpRadii, ...
    'resharpRegularization', cfg.resharpRegularization);
stepTimer = tic;
fmriPhase.util.log(cfg, 'Computing harmonic off-resonance fields. This includes RESHARP and may take a while.');
fields = fmriPhase.phase.computeHarmonicFields(data.phase, mask0, fieldParams);
fmriPhase.util.logElapsed(cfg, 'Harmonic field computation', stepTimer);

switch cfg.inputImage
    case "harmonicField"
        svdInput = fields.harmonicField;
    case "totalField"
        svdInput = fields.totalField;
    case "localField"
        svdInput = fields.localField;
    otherwise
        error('fmriPhase:pipeline:UnsupportedInputImage', ...
            'Unsupported inputImage "%s".', cfg.inputImage);
end

stepTimer = tic;
fmriPhase.util.log(cfg, 'Extracting respiratory SVD component.');
respiration = fmriPhase.respiration.extractSvdComponent( ...
    svdInput, data.TR, ...
    'interpolationFactor', cfg.interpolationFactor, ...
    'nComponents', cfg.nSvdComponents);
fmriPhase.util.logElapsed(cfg, 'SVD respiratory extraction', stepTimer);

stepTimer = tic;
fmriPhase.util.log(cfg, 'Projecting respiratory field onto harmonic basis.');
coefficients = fmriPhase.harmonics.projectTimecourse( ...
    respiration.respiratoryVolume, fields.mask, data.voxelSize, cfg.harmonicOrder);
regressors = normalize([coefficients, respiration.componentTimecourses], 1);
regressorNames = fmriPhase.regressors.defaultNames( ...
    size(coefficients, 2), size(respiration.componentTimecourses, 2));
fmriPhase.util.logElapsed(cfg, 'Harmonic projection', stepTimer);

result = struct();
result.config = cfg;
result.data = data;
result.mask = fields.mask;
result.fields = fields;
result.respiration = respiration;
result.harmonicCoefficients = coefficients;
result.regressors = regressors;
result.regressorNames = regressorNames;
result.TR = data.TR;
result.runtimeSeconds = toc(pipelineTimer);
fmriPhase.util.logElapsed(cfg, 'Pipeline', pipelineTimer);

if cfg.saveRegressors
    stepTimer = tic;
    fmriPhase.util.log(cfg, 'Saving MATLAB regressor archive.');
    reg = regressors;
    names = regressorNames;
    harmonicCoefficients = coefficients;
    respirationSummary = rmfield(respiration, 'respiratoryVolume');
    save(fullfile(cfg.outputDir, cfg.regressorFilename), ...
        'reg', 'names', 'harmonicCoefficients', 'respirationSummary', 'cfg');
    fmriPhase.util.logElapsed(cfg, 'MATLAB archive save', stepTimer);
end

end
