function outputs = generateReport(result, regressorOutputs, cfg)
%GENERATEREPORT Write QC figures and a text summary for one pipeline run.

arguments
    result struct
    regressorOutputs struct
    cfg struct
end

qcDir = string(cfg.qcDir);
if strlength(qcDir) == 0
    qcDir = fullfile(string(cfg.outputDir), string(cfg.regressorPrefix) + "_qc");
end

if ~isfolder(qcDir)
    mkdir(qcDir);
end

outputs = struct();
outputs.directory = qcDir;

summary = buildSummary(result, cfg);
outputs.summary = fullfile(qcDir, "qc_summary.txt");
writeSummary(outputs.summary, summary);

outputs.regressors = fullfile(qcDir, "regressors.png");
plotRegressors(result, outputs.regressors);

outputs.svdSelection = fullfile(qcDir, "svd_selection.png");
plotSvdSelection(result, outputs.svdSelection);

outputs.respiratoryPhase = fullfile(qcDir, "respiratory_phase.png");
plotRespiratoryPhase(result, outputs.respiratoryPhase);

outputs.breathingField = fullfile(qcDir, "breathing_off_resonance_field.png");
plotBreathingField(result, outputs.breathingField);

outputs.regressorCorrelation = fullfile(qcDir, "regressor_correlation.png");
plotRegressorCorrelation(result, outputs.regressorCorrelation);

if strlength(string(cfg.motionFile)) > 0
    motionQc = fmriPhase.qc.motionQc(result, cfg.motionFile, cfg.motionModel);
    outputs.motionCorrelation = fullfile(qcDir, "motion_correlation.png");
    plotMotionCorrelation(motionQc, outputs.motionCorrelation);
    outputs.motionOverlay = fullfile(qcDir, "selected_svd_vs_motion.png");
    plotMotionOverlay(result, motionQc, outputs.motionOverlay);
    summary.MotionFile = char(cfg.motionFile);
    summary.SelectedSvdMotionR2 = motionQc.selectedSvdMotionR2;
    writeSummary(outputs.summary, summary);
end

outputs.regressorOutputs = regressorOutputs;

end

function summary = buildSummary(result, cfg)
summary = struct();
summary.PhaseFile = char(cfg.phaseFile);
summary.MagnitudeFile = char(cfg.magnitudeFile);
summary.OutputDir = char(cfg.outputDir);
summary.TR = result.TR;
summary.NumberOfRegressors = size(result.regressors, 2);
summary.NumberOfVolumes = size(result.regressors, 1);
summary.SelectedSvdComponent = result.respiration.selectedComponent;
summary.RuntimeSeconds = getOptionalField(result, 'runtimeSeconds', NaN);
summary.InputImage = char(cfg.inputImage);
summary.HarmonicOrder = cfg.harmonicOrder;
end

function writeSummary(file, summary)
fid = fopen(file, 'w');
if fid < 0
    error('fmriPhase:qc:OpenFile', 'Could not open "%s" for writing.', file);
end
cleanup = onCleanup(@() fclose(fid));
fields = fieldnames(summary);
for idx = 1:numel(fields)
    value = summary.(fields{idx});
    if isnumeric(value)
        valueText = mat2str(value);
    else
        valueText = char(string(value));
    end
    fprintf(fid, '%s: %s\n', fields{idx}, valueText);
end
end

function plotRegressors(result, file)
fig = newHiddenFigure();
plot(result.regressors, 'LineWidth', 0.8);
title('fMRI phase-derived regressors');
xlabel('Volume');
ylabel('Normalized amplitude');
legend(cellstr(result.regressorNames), 'Interpreter', 'none', 'Location', 'eastoutside');
saveAndClose(fig, file);
end

function plotSvdSelection(result, file)
fig = newHiddenFigure();
bar(result.respiration.selectionScores);
hold on;
xline(result.respiration.selectedComponent, 'r', 'LineWidth', 2);
title('SVD respiratory component selection scores');
xlabel('SVD component');
ylabel('L1 norm of temporal derivative');
saveAndClose(fig, file);
end

function plotRespiratoryPhase(result, file)
fig = newHiddenFigure();
tiledlayout(3, 1);
nexttile;
plot(result.respiration.selectedTimecourse, 'LineWidth', 1.2);
title('Selected SVD respiratory timecourse');
xlabel('Volume');
ylabel('Amplitude');
nexttile;
plot(result.respiration.respiratoryPhase, 'LineWidth', 1.2);
ylim([0 2*pi]);
title('Respiratory phase');
xlabel('Volume');
ylabel('Phase [rad]');
nexttile;
plot(result.respiration.respiratoryPhaseSin, 'LineWidth', 1.2);
hold on;
plot(result.respiration.respiratoryPhaseCos, 'LineWidth', 1.2);
legend({'sin(phase)', 'cos(phase)'});
title('Circular respiratory phase terms');
xlabel('Volume');
ylabel('Amplitude');
saveAndClose(fig, file);
end

function plotBreathingField(result, file)
respiration = result.respiration;
[~, lowIdx] = min(respiration.selectedTimecourse);
[~, highIdx] = max(respiration.selectedTimecourse);
lowField = respiration.respiratoryVolume(:, :, :, lowIdx);
highField = respiration.respiratoryVolume(:, :, :, highIdx);
peakToPeakField = highField - lowField;

sliceStd = squeeze(std(peakToPeakField, 0, [1 2], 'omitnan'));
[~, fieldSlice] = max(sliceStd);
fieldLimit = max(abs(peakToPeakField(:, :, fieldSlice)), [], 'all');
if fieldLimit == 0 || ~isfinite(fieldLimit)
    fieldLimit = 1;
end

fig = newHiddenFigure();
tiledlayout(1, 3);
nexttile;
imagesc(lowField(:, :, fieldSlice));
axis image off;
colorbar;
title(sprintf('Low field, slice %d', fieldSlice));
nexttile;
imagesc(highField(:, :, fieldSlice));
axis image off;
colorbar;
title(sprintf('High field, slice %d', fieldSlice));
nexttile;
imagesc(peakToPeakField(:, :, fieldSlice));
axis image off;
colorbar;
clim(gca, [-fieldLimit fieldLimit]);
title('High - low field');
colormap(gca, turbo);
sgtitle('Breathing off-resonance field diagnostic');
saveAndClose(fig, file);
end

function plotRegressorCorrelation(result, file)
fig = newHiddenFigure();
imagesc(corr(result.regressors, 'Rows', 'pairwise'));
axis image;
colorbar;
clim(gca, [-1 1]);
title('Regressor correlation matrix');
xticks(1:numel(result.regressorNames));
yticks(1:numel(result.regressorNames));
xticklabels(result.regressorNames);
yticklabels(result.regressorNames);
xtickangle(90);
set(gca, 'TickLabelInterpreter', 'none');
saveAndClose(fig, file);
end

function plotMotionCorrelation(motionQc, file)
fig = newHiddenFigure();
imagesc(motionQc.correlation);
colorbar;
clim(gca, [-1 1]);
title('Phase regressors vs motion model');
xlabel('Motion model column');
ylabel('Phase regressor');
yticks(1:numel(motionQc.regressorNames));
yticklabels(motionQc.regressorNames);
set(gca, 'TickLabelInterpreter', 'none');
saveAndClose(fig, file);
end

function plotMotionOverlay(result, motionQc, file)
fig = newHiddenFigure();
selected = normalize(result.respiration.selectedTimecourse, 1);
yyaxis left;
plot(selected, 'LineWidth', 1.2);
ylabel('Selected SVD timecourse');
yyaxis right;
plot(motionQc.framewiseDisplacementApprox, 'LineWidth', 1.0);
ylabel('Approximate framewise displacement');
title(sprintf('Selected SVD vs motion, R^2 = %.3f', motionQc.selectedSvdMotionR2));
xlabel('Volume');
saveAndClose(fig, file);
end

function fig = newHiddenFigure()
fig = figure('Visible', 'off', 'Color', 'w');
set(fig, 'Position', [100 100 1200 700]);
end

function saveAndClose(fig, file)
exportgraphics(fig, file, 'Resolution', 150);
close(fig);
end

function value = getOptionalField(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end
