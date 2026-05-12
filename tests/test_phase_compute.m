classdef test_phase_compute < matlab.unittest.TestCase
    properties
        MockDir
    end

    methods (TestMethodSetup)
        function addMockQsmFunctions(testCase)
            testCase.MockDir = tempname;
            mkdir(testCase.MockDir);
            testCase.writeMockFunction('unwrapLaplacian.m', [
                "function uphas = unwrapLaplacian(phas, mask, voxelSize)" newline ...
                "% Minimal test double for QSM.m unwrapLaplacian." newline ...
                "uphas = phas + double(mask) * 0 + sum(voxelSize) * 0;" newline ...
                "end" newline]);
            testCase.writeMockFunction('resharp.m', [
                "function [localField, maskOut] = resharp(totalField, mask, voxelSize, radii, lambda)" newline ...
                "% Minimal test double for QSM.m resharp." newline ...
                "localField = totalField * 0.25 + sum(voxelSize) * 0 + sum(radii) * 0 + lambda * 0;" newline ...
                "maskOut = mask;" newline ...
                "end" newline]);
            addpath(testCase.MockDir, '-begin');
        end
    end

    methods (TestMethodTeardown)
        function removeMockQsmFunctions(testCase)
            rmpath(testCase.MockDir);
            rmdir(testCase.MockDir, 's');
        end
    end

    methods (Test)
        function acceptsPipelineParameterStruct(testCase)
            phaseData = int16(zeros(2, 2, 2, 3));
            mask = true(2, 2, 2);
            params = struct( ...
                'echoTime', 0.030, ...
                'fieldStrength', 3.0, ...
                'gyromagneticRatio', 267.513, ...
                'voxelSize', [2 2 2], ...
                'resharpRadii', [], ...
                'resharpRegularization', 0.05);

            fields = fmriPhase.phase.computeHarmonicFields(phaseData, mask, params);

            testCase.verifyTrue(isfield(fields, 'phaseRadians'));
            testCase.verifyTrue(isfield(fields, 'unwrappedPhase'));
            testCase.verifyTrue(isfield(fields, 'totalField'));
            testCase.verifyTrue(isfield(fields, 'localField'));
            testCase.verifyTrue(isfield(fields, 'harmonicField'));
            testCase.verifyEqual(fields.mask, mask);
            testCase.verifySize(fields.harmonicField, size(phaseData));
        end
    end

    methods (Access = private)
        function writeMockFunction(testCase, filename, contents)
            fid = fopen(fullfile(testCase.MockDir, filename), 'w');
            testCase.assertGreaterThan(fid, 0);
            cleanup = onCleanup(@() fclose(fid));
            fprintf(fid, '%s', contents);
        end
    end
end
