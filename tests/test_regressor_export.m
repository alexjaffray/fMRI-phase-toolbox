classdef test_regressor_export < matlab.unittest.TestCase
    methods (Test)
        function writesPackageCompatibleFiles(testCase)
            outDir = tempname;
            mkdir(outDir);
            cleanup = onCleanup(@() rmdir(outDir, 's'));

            result = struct();
            result.regressors = reshape(1:12, [4, 3]);
            result.regressorNames = ["phase_harmonic_01", "phase_harmonic_02", "phase_svd_01"];
            result.TR = 1.15;
            result.config = fmriPhase.config( ...
                'phaseFile', "phase.nii", ...
                'magnitudeFile', "mag.nii");
            result.respiration = struct();
            result.respiration.selectedComponent = 1;
            result.respiration.selectionScores = [1; 2; 3];
            result.respiration.singularValues = [10; 5; 1];
            result.respiration.selectedTimecourse = [0; 1; -1; 2];
            result.respiration.respiratoryPhase = [0; 1; 2; 3];
            result.respiration.respiratoryPhaseSin = sin(result.respiration.respiratoryPhase);
            result.respiration.respiratoryPhaseCos = cos(result.respiration.respiratoryPhase);
            result.respiration.respiratoryVolume = reshape(1:32, [2, 2, 2, 4]);

            outputs = fmriPhase.regressors.export(result, outDir, "sub-01_task-rest");

            testCase.verifyTrue(isfile(outputs.matlab));
            testCase.verifyTrue(isfile(outputs.bids.tsv));
            testCase.verifyTrue(isfile(outputs.bids.json));
            testCase.verifyTrue(isfile(outputs.spm.mat));
            testCase.verifyTrue(isfile(outputs.spm.txt));
            testCase.verifyTrue(isfile(outputs.afni));
            testCase.verifyTrue(isfile(outputs.conn.tsv));
            testCase.verifyTrue(isfile(outputs.conn.mat));
            testCase.verifyTrue(isfolder(outputs.conn.perCovariateDir));

            spm = load(outputs.spm.mat);
            testCase.verifyEqual(spm.R, result.regressors);

            exported = load(outputs.matlab);
            testCase.verifyTrue(isfield(exported, 'svdDiagnostics'));
            testCase.verifyTrue(isfield(exported.svdDiagnostics, 'RespiratoryPhase'));
            testCase.verifyTrue(isfield(exported.svdDiagnostics, 'RespiratoryPhaseSin'));
            testCase.verifyTrue(isfield(exported.svdDiagnostics, 'RespiratoryPhaseCos'));
            testCase.verifyTrue(isfield(exported.svdDiagnostics, 'BreathingFieldPeakToPeak'));
            testCase.verifySize(exported.svdDiagnostics.BreathingFieldPeakToPeak, [2, 2, 2]);

            tsvText = fileread(outputs.bids.tsv);
            testCase.verifyTrue(startsWith(tsvText, 'phase_harmonic_01'));
        end
    end
end
