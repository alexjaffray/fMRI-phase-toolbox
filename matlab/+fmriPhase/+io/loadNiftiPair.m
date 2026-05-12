function data = loadNiftiPair(phaseFile, magnitudeFile, rotate90)
%LOADNIFTIPAIR Load phase and magnitude fMRI NIfTI files.

arguments
    phaseFile {mustBeTextScalar}
    magnitudeFile {mustBeTextScalar}
    rotate90 (1,1) logical = false
end

data = struct();
data.phaseFile = string(phaseFile);
data.magnitudeFile = string(magnitudeFile);
data.phase = niftiread(data.phaseFile);
data.magnitude = niftiread(data.magnitudeFile);
data.info = niftiinfo(data.phaseFile);

if rotate90
    data.phase = rot90(data.phase);
    data.magnitude = rot90(data.magnitude);
end

if ndims(data.phase) ~= 4
    error('fmriPhase:io:PhaseDimensions', ...
        'Phase image must be 4-D, got %d dimensions.', ndims(data.phase));
end

if ~isequal(size(data.phase), size(data.magnitude))
    error('fmriPhase:io:ImageSizeMismatch', ...
        'Phase and magnitude images must have identical dimensions.');
end

data.imageSize = size(data.phase);
data.voxelSize = data.info.PixelDimensions(1:3);
data.TR = data.info.PixelDimensions(4);

end
