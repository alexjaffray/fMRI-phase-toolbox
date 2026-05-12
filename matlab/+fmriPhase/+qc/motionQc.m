function qc = motionQc(result, motionFile, motionModel)
%MOTIONQC Compare phase regressors and selected SVD component with motion.

arguments
    result struct
    motionFile {mustBeTextScalar}
    motionModel (1,1) string = "motion+derivatives"
end

motion = readmatrix(motionFile);
if size(motion, 1) ~= size(result.regressors, 1)
    error('fmriPhase:qc:MotionLengthMismatch', ...
        'Motion file has %d rows, but regressors have %d rows.', ...
        size(motion, 1), size(result.regressors, 1));
end

switch lower(motionModel)
    case "motion"
        Xmotion = motion;
    case "motion+derivatives"
        Xmotion = [motion, [zeros(1, size(motion, 2)); diff(motion)]];
    case "motion+derivatives+squares"
        derivatives = [zeros(1, size(motion, 2)); diff(motion)];
        Xmotion = [motion, derivatives, motion.^2, derivatives.^2];
    otherwise
        error('fmriPhase:qc:UnsupportedMotionModel', ...
            'Unsupported motionModel "%s".', motionModel);
end

Xmotion = normalize(Xmotion, 1);
selected = normalize(result.respiration.selectedTimecourse, 1);
design = [ones(size(Xmotion, 1), 1), Xmotion];
fit = design \ selected;
prediction = design * fit;
residual = selected - prediction;
denominator = sum((selected - mean(selected, 'omitnan')).^2, 'omitnan');
if denominator == 0
    selectedR2 = NaN;
else
    selectedR2 = 1 - sum(residual.^2, 'omitnan') / denominator;
end

qc = struct();
qc.motion = motion;
qc.motionModel = Xmotion;
qc.correlation = corr(result.regressors, Xmotion, 'Rows', 'pairwise');
qc.regressorNames = result.regressorNames;
qc.selectedSvdMotionR2 = selectedR2;
qc.framewiseDisplacementApprox = sum(abs([zeros(1, size(motion, 2)); diff(motion)]), 2);

end
