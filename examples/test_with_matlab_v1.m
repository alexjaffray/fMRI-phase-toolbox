% Example MATLAB script to test loading and visualizing fMRI phase-derived regressors including my own testing paths. 
% Edit these to make them work for you

close all;
clc;

cd('/srv/data/ajaffray/fMRI-phase-toolbox-v2')
addpath('/srv/data/ajaffray/fMRI-phase-toolbox-v2/matlab')
addpath('/srv/data/ajaffray/fMRI-phase-toolbox-v2/matlab/vendor/tapas_physio')

cfg = fmriPhase.config();

cfg.phaseFile = "/srv/data/ajaffray/fMRI-phase-toolbox/template/sub-M02_ses-2122pre_task-rest_run-1_part-phase_bold.nii";
cfg.magnitudeFile = "/srv/data/ajaffray/fMRI-phase-toolbox/template/sub-M02_ses-2122pre_task-rest_run-1_part-mag_bold.nii";
cfg.outputDir = "/srv/data/ajaffray/fMRI-phase-toolbox-v2/outputs/sub-M02_ses-2122pre_run-1";
cfg.regressorPrefix = "sub-M02_ses-2122pre_task-rest_run-1";

cfg.qsmSetup = "/srv/data/ajaffray/QSM/addpathqsm.m";

cfg.verbose = true;
cfg.generateQc = true;

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
cfg.exportFormats = ["matlab", "bids", "spm", "afni", "conn"];

outputs = fmriPhase.pipeline.generateRegressors(cfg);

disp(outputs)

cfg.outputDir

S = load(outputs.matlab);
size(S.regressors)
disp(S.regressorNames')

svdIdx = startsWith(string(S.regressorNames), "phase_svd_");
svdRegressors = S.regressors(:, svdIdx);
svdNames = string(S.regressorNames(svdIdx));

if isfield(S, 'svdDiagnostics')
    fprintf('Selected SVD component: %d\n', S.svdDiagnostics.SelectedComponent);
    disp('SVD component selection scores:');
    disp(S.svdDiagnostics.SelectionScores(:)');
    if isfield(S.svdDiagnostics, 'RespiratoryPhase')
        disp('Respiratory phase range [rad]:');
        disp([min(S.svdDiagnostics.RespiratoryPhase), max(S.svdDiagnostics.RespiratoryPhase)]);
    else
        warning(['This MAT file has SVD diagnostics but no respiratory phase. ', ...
            'Clear functions and rerun fmriPhase.pipeline.generateRegressors(cfg) ', ...
            'to regenerate the outputs with the current v2 exporter.']);
    end
end
%%
figure;
plot(S.regressors);
legend(S.regressorNames, 'Interpreter', 'none');
title('fMRI phase-derived regressors');

figure;
plot(svdRegressors, 'LineWidth', 1.2);
legend(svdNames, 'Interpreter', 'none');
title('SVD regressors');
xlabel('Volume');
ylabel('Normalized amplitude');

if isfield(S, 'svdDiagnostics')
    figure;
    bar(S.svdDiagnostics.SelectionScores);
    hold on;
    xline(S.svdDiagnostics.SelectedComponent, 'r', 'LineWidth', 2);
    title('SVD respiratory component selection scores');
    xlabel('SVD component');
    ylabel('L1 norm of temporal derivative');
end

if isfield(S, 'svdDiagnostics') && isfield(S.svdDiagnostics, 'RespiratoryPhase')
    figure;
    subplot(3,1,1);
    plot(S.svdDiagnostics.SelectedTimecourse, 'LineWidth', 1.2);
    title('Selected SVD respiratory timecourse');
    xlabel('Volume');
    ylabel('Amplitude');

    subplot(3,1,2);
    plot(S.svdDiagnostics.RespiratoryPhase, 'LineWidth', 1.2);
    title('Respiratory phase from selected SVD component');
    xlabel('Volume');
    ylabel('Phase [rad]');
    ylim([0 2*pi]);

    subplot(3,1,3);
    plot(S.svdDiagnostics.RespiratoryPhaseSin, 'LineWidth', 1.2);
    hold on;
    plot(S.svdDiagnostics.RespiratoryPhaseCos, 'LineWidth', 1.2);
    legend({'sin(phase)', 'cos(phase)'});
    title('Circular respiratory phase terms');
    xlabel('Volume');
    ylabel('Amplitude');
end

if isfield(S, 'svdDiagnostics') && isfield(S.svdDiagnostics, 'BreathingFieldPeakToPeak') ...
        && ~isempty(S.svdDiagnostics.BreathingFieldPeakToPeak)
    peakToPeakField = S.svdDiagnostics.BreathingFieldPeakToPeak;
    lowField = S.svdDiagnostics.BreathingFieldLow;
    highField = S.svdDiagnostics.BreathingFieldHigh;

    fieldStdBySlice = squeeze(std(peakToPeakField, 0, [1 2], 'omitnan'));
    [~, fieldSlice] = max(fieldStdBySlice);

    fieldLimit = max(abs(peakToPeakField(:, :, fieldSlice)), [], 'all');
    if fieldLimit == 0 || ~isfinite(fieldLimit)
        fieldLimit = 1;
    end

    figure;
    tiledlayout(1, 3);

    nexttile;
    imagesc(lowField(:, :, fieldSlice));
    axis image off;
    colorbar;
    title(sprintf('Low breathing field, slice %d', fieldSlice));

    nexttile;
    imagesc(highField(:, :, fieldSlice));
    axis image off;
    colorbar;
    title(sprintf('High breathing field, slice %d', fieldSlice));

    nexttile;
    imagesc(peakToPeakField(:, :, fieldSlice));
    axis image off;
    colorbar;
    clim(gca, [-fieldLimit fieldLimit]);
    title('Breathing off-resonance field');
    colormap(gca, turbo);

    sgtitle(sprintf('Selected SVD component %d: high - low field', ...
        S.svdDiagnostics.SelectedComponent));
else
    warning(['No breathing field maps found in this MAT file. ', ...
        'Regenerate outputs with the current v2 exporter to visualize the off-resonance field.']);
end
