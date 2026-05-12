function outputs = generateRegressors(cfg)
%GENERATEREGRESSORS Run the pipeline and export shareable nuisance regressors.
%
% outputs = fmriPhase.pipeline.generateRegressors(cfg)
%
% This is the collaborator-facing entry point. It produces one regressor
% matrix and writes files that can be imported by common fMRI packages.

cfg = fmriPhase.validateConfig(cfg);
if strlength(string(cfg.regressorPrefix)) == 0
    [~, baseName] = fileparts(string(cfg.phaseFile));
    cfg.regressorPrefix = regexprep(baseName, '_part-phase_bold$', '');
end

result = fmriPhase.pipeline.run(cfg);
fmriPhase.util.log(cfg, 'Exporting package-compatible regressor files.');
outputs = fmriPhase.regressors.export(result, cfg.outputDir, cfg.regressorPrefix, ...
    'formats', cfg.exportFormats);

if cfg.generateQc
    fmriPhase.util.log(cfg, 'Generating QC report.');
    outputs.qc = fmriPhase.qc.generateReport(result, outputs, cfg);
end

end
