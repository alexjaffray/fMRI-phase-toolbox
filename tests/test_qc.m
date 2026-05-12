classdef test_qc < matlab.unittest.TestCase
    methods (Test)
        function motionQcComputesSelectedSvdR2(testCase)
            outDir = tempname;
            mkdir(outDir);
            cleanup = onCleanup(@() rmdir(outDir, 's'));

            motionFile = fullfile(outDir, 'motion.txt');
            motion = [(1:8).', sin((1:8).'), cos((1:8).'), zeros(8, 3)];
            writematrix(motion, motionFile, 'Delimiter', 'tab');

            result = minimalResult(8);
            qc = fmriPhase.qc.motionQc(result, motionFile, "motion+derivatives");

            testCase.verifySize(qc.correlation, [size(result.regressors, 2), 12]);
            testCase.verifySize(qc.framewiseDisplacementApprox, [8, 1]);
            testCase.verifyTrue(isfinite(qc.selectedSvdMotionR2));
        end

        function reportWritesExpectedFiles(testCase)
            outDir = tempname;
            mkdir(outDir);
            cleanup = onCleanup(@() rmdir(outDir, 's'));

            result = minimalResult(8);
            cfg = fmriPhase.config( ...
                'phaseFile', "phase.nii", ...
                'magnitudeFile', "mag.nii", ...
                'outputDir', outDir, ...
                'regressorPrefix', "sub-01_task-rest", ...
                'generateQc', true, ...
                'verbose', false);

            outputs = fmriPhase.qc.generateReport(result, struct(), cfg);

            testCase.verifyTrue(isfile(outputs.summary));
            testCase.verifyTrue(isfile(outputs.regressors));
            testCase.verifyTrue(isfile(outputs.svdSelection));
            testCase.verifyTrue(isfile(outputs.respiratoryPhase));
            testCase.verifyTrue(isfile(outputs.breathingField));
            testCase.verifyTrue(isfile(outputs.regressorCorrelation));
        end
    end
end

function result = minimalResult(nTimepoints)
names = fmriPhase.regressors.defaultNames(2, 2);
time = (0:(nTimepoints - 1)).';
selected = sin(2 * pi * time / nTimepoints);

result = struct();
result.config = fmriPhase.config('verbose', false);
result.TR = 1.0;
result.runtimeSeconds = 1.23;
result.regressors = normalize([selected, cos(2 * pi * time / nTimepoints), ...
    selected.^2, cos(2 * pi * time / nTimepoints).^2], 1);
result.regressorNames = names;

respiration = struct();
respiration.selectedComponent = 1;
respiration.selectionScores = [2; 1];
respiration.singularValues = [10; 5];
respiration.selectedTimecourse = selected;
respiration.respiratoryPhase = mod(atan2(sin(time), cos(time)), 2 * pi);
respiration.respiratoryPhaseSin = sin(respiration.respiratoryPhase);
respiration.respiratoryPhaseCos = cos(respiration.respiratoryPhase);
respiration.respiratoryVolume = zeros(4, 4, 3, nTimepoints);
for idx = 1:nTimepoints
    respiration.respiratoryVolume(:, :, :, idx) = idx;
end
result.respiration = respiration;
end
