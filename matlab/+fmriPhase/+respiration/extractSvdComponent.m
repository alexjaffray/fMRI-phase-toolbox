function result = extractSvdComponent(inputImage, TR, opts)
%EXTRACTSVDCOMPONENT Extract respiratory-dominant image component by SVD.

arguments
    inputImage
    TR (1,1) double {mustBePositive}
    opts.interpolationFactor (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.nComponents (1,1) double {mustBeInteger, mustBePositive} = 5
end

imageSize = size(inputImage);
if numel(imageSize) ~= 4
    error('fmriPhase:respiration:InputDimensions', ...
        'inputImage must be 4-D.');
end

nTimepoints = imageSize(4);
nComponents = min([opts.nComponents, nTimepoints, prod(imageSize(1:3))]);
matrix = reshape(inputImage, prod(imageSize(1:3)), nTimepoints);
[U, S, V] = svd(double(matrix), "econ");

scores = zeros(nComponents, 1);
componentVolumes = cell(nComponents, 1);
componentTimecourses = zeros(nTimepoints, nComponents);

for componentIdx = 1:nComponents
    [componentVolumes{componentIdx}, raw] = fmriPhase.respiration.recomposeSvd( ...
        U, S, V, componentIdx, imageSize(1:3));
    componentTimecourses(:, componentIdx) = raw(componentIdx, :).';
    scores(componentIdx) = norm(diff(componentTimecourses(:, componentIdx)), 1);
end

[~, selectedComponent] = max(scores);
selectedVolume = componentVolumes{selectedComponent};

if opts.interpolationFactor > 1
    selectedVolume = interpft(selectedVolume, nTimepoints * opts.interpolationFactor, 4);
    componentTimecoursesOut = interpft(componentTimecourses, ...
        nTimepoints * opts.interpolationFactor, 1);
else
    componentTimecoursesOut = componentTimecourses;
end

time = (0:(size(selectedVolume, 4) - 1)).' * TR / opts.interpolationFactor;
selectedTimecourse = componentTimecoursesOut(:, selectedComponent);
respiratoryPhase = fmriPhase.respiration.estimatePhase(selectedTimecourse, time);

result = struct();
result.respiratoryVolume = selectedVolume;
result.time = time;
result.zeroOrderVolume = componentVolumes{1};
result.componentTimecourses = componentTimecoursesOut;
result.rawComponentTimecourses = componentTimecourses;
result.selectionScores = scores;
result.selectedComponent = selectedComponent;
result.selectedTimecourse = selectedTimecourse;
result.respiratoryPhase = respiratoryPhase;
result.respiratoryPhaseSin = sin(respiratoryPhase);
result.respiratoryPhaseCos = cos(respiratoryPhase);
result.singularValues = diag(S);

end
