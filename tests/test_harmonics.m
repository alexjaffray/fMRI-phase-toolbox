classdef test_harmonics < matlab.unittest.TestCase
    methods (Test)
        function designMatrixHasExpectedOrder3Columns(testCase)
            mask = true(4, 5, 3);
            basis = fmriPhase.harmonics.designMatrix(mask, [1 1 1], 3);
            testCase.verifySize(basis, [numel(mask), 16]);
        end

        function projectionRecoversConstantTerm(testCase)
            mask = true(4, 5, 3);
            image = ones(4, 5, 3, 2);
            coeffs = fmriPhase.harmonics.projectTimecourse(image, mask, [1 1 1], 1);
            testCase.verifyEqual(coeffs(:, 1), [1; 1], 'AbsTol', 1e-10);
        end
    end
end
