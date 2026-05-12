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

fields = struct();
fields.phaseRadians = double(phaseData) ./ 2048 * pi - pi;
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
