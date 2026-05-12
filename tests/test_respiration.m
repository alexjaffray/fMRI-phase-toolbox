classdef test_respiration < matlab.unittest.TestCase
    methods (Test)
        function extractsDominantChangingComponent(testCase)
            t = linspace(0, 2 * pi, 12);
            image = zeros(3, 3, 2, numel(t));
            image(1, 1, 1, :) = sin(t);
            image(2, 2, 1, :) = 0.1 * cos(t);

            result = fmriPhase.respiration.extractSvdComponent(image, 1.0, ...
                'nComponents', 2, 'interpolationFactor', 1);

            testCase.verifyGreaterThanOrEqual(result.selectedComponent, 1);
            testCase.verifyLessThanOrEqual(result.selectedComponent, 2);
            testCase.verifySize(result.time, [12, 1]);
            testCase.verifySize(result.componentTimecourses, [12, 2]);
            testCase.verifySize(result.selectedTimecourse, [12, 1]);
            testCase.verifySize(result.respiratoryPhase, [12, 1]);
            testCase.verifyGreaterThanOrEqual(min(result.respiratoryPhase), 0);
            testCase.verifyLessThan(max(result.respiratoryPhase), 2 * pi + eps);
            testCase.verifySize(result.respiratoryPhaseSin, [12, 1]);
            testCase.verifySize(result.respiratoryPhaseCos, [12, 1]);
        end

        function interpolationKeepsTimecoursesAligned(testCase)
            t = linspace(0, 2 * pi, 12);
            image = zeros(3, 3, 2, numel(t));
            image(1, 1, 1, :) = sin(t);

            result = fmriPhase.respiration.extractSvdComponent(image, 1.0, ...
                'nComponents', 1, 'interpolationFactor', 2);

            testCase.verifySize(result.time, [24, 1]);
            testCase.verifySize(result.componentTimecourses, [24, 1]);
            testCase.verifySize(result.respiratoryPhase, [24, 1]);
        end
    end
end
