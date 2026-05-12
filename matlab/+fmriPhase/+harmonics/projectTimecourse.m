function coefficients = projectTimecourse(image4d, mask, voxelSize, order)
%PROJECTTIMECOURSE Project 4-D field data onto harmonic basis terms.

arguments
    image4d
    mask
    voxelSize (1,3) double {mustBePositive} = [1 1 1]
    order (1,1) double {mustBeInteger, mustBePositive} = 3
end

if ndims(image4d) ~= 4
    error('fmriPhase:harmonics:InputDimensions', ...
        'image4d must be 4-D.');
end

basis = fmriPhase.harmonics.designMatrix(mask, voxelSize, order);
nTimepoints = size(image4d, 4);
coefficients = zeros(nTimepoints, size(basis, 2));

for timeIdx = 1:nTimepoints
    volume = image4d(:, :, :, timeIdx);
    coefficients(timeIdx, :) = basis \ volume(mask);
end

end
