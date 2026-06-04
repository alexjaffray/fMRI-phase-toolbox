function fields = computeHarmonicFields(phaseData, mask, params)
%COMPUTEHARMONICFIELDS Unwrap phase and isolate harmonic background fields.

arguments
    phaseData
    mask
    params (1,1) struct
end

requiredFields = ["echoTime", "fieldStrength", "gyromagneticRatio", "voxelSize"];
for fieldName = requiredFields
    if ~isfield(params, fieldName)
        error('fmriPhase:phase:MissingParameter', ...
            'Missing phase parameter "%s".', fieldName);
    end
end

if ~isfield(params, 'resharpRadii')
    params.resharpRadii = [];
end

if ~isfield(params, 'resharpRegularization')
    params.resharpRegularization = 0.05;
end

if ~isfield(params, 'phaseInputUnits')
    params.phaseInputUnits = "scaled";
end

if ~isfield(params, 'phaseScale')
    params.phaseScale = pi / 2048;
end

if ~isfield(params, 'phaseOffset')
    params.phaseOffset = -pi;
end

if ~isfield(params, 'phaseRangeTolerance')
    params.phaseRangeTolerance = 0.1;
end

fields = struct();
[fields.phaseRadians, fields.phaseScaling] = convertPhaseInput(phaseData, params);
fields.unwrappedPhase = unwrapLaplacian(fields.phaseRadians, mask, params.voxelSize);

scale = params.fieldStrength * params.gyromagneticRatio * params.echoTime;
fields.totalField = fields.unwrappedPhase ./ scale;

if isempty(params.resharpRadii)
    maxVoxel = max(params.voxelSize);
    params.resharpRadii = 9:-2 * maxVoxel:2 * maxVoxel;
end

[fields.localField, fields.mask] = resharp( ...
    fields.totalField, mask, params.voxelSize, ...
    params.resharpRadii, params.resharpRegularization);
fields.harmonicField = fields.totalField - fields.localField;

end

function [phaseRadians, scaling] = convertPhaseInput(phaseData, params)
phaseValues = double(phaseData);
finiteValues = phaseValues(isfinite(phaseValues));
if isempty(finiteValues)
    phaseRange = NaN;
else
    phaseRange = max(finiteValues(:)) - min(finiteValues(:));
end

requestedUnits = string(params.phaseInputUnits);
detectedUnits = requestedUnits;
if requestedUnits == "auto"
    if isfinite(phaseRange) && abs(phaseRange - 2 * pi) <= params.phaseRangeTolerance
        detectedUnits = "radians";
    else
        detectedUnits = "scaled";
    end
end

switch detectedUnits
    case "radians"
        phaseRadians = phaseValues;
        scale = 1;
        offset = 0;
    case "scaled"
        if isempty(params.phaseScale)
            error('fmriPhase:phase:PhaseScaleRequired', ...
                ['phaseScale must be specified when phaseInputUnits is "scaled" ', ...
                'or when auto-detection does not find a 2*pi phase range.']);
        end
        scale = params.phaseScale;
        offset = params.phaseOffset;
        phaseRadians = phaseValues .* scale + offset;
    otherwise
        error('fmriPhase:phase:UnsupportedPhaseInputUnits', ...
            'Unsupported phaseInputUnits "%s".', requestedUnits);
end

scaling = struct( ...
    'requestedUnits', requestedUnits, ...
    'detectedUnits', detectedUnits, ...
    'inputRange', phaseRange, ...
    'rangeTolerance', params.phaseRangeTolerance, ...
    'scale', scale, ...
    'offset', offset);
end
