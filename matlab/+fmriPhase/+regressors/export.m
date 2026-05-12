function outputs = export(result, outputDir, prefix, opts)
%EXPORT Write regressor outputs for MATLAB, BIDS-style, SPM, AFNI, and CONN.

arguments
    result struct
    outputDir {mustBeTextScalar}
    prefix {mustBeTextScalar}
    opts.formats = ["matlab", "bids", "spm", "afni", "conn"]
end

outputDir = string(outputDir);
prefix = string(prefix);
formats = lower(string(opts.formats));

if ~isfolder(outputDir)
    mkdir(outputDir);
end

R = result.regressors;
names = string(result.regressorNames);
if size(R, 2) ~= numel(names)
    error('fmriPhase:regressors:NameCountMismatch', ...
        'Number of regressor names must match number of matrix columns.');
end

outputs = struct();

if any(formats == "matlab")
    outputs.matlab = writeMatlab(outputDir, prefix, result, R, names);
end

if any(formats == "bids")
    outputs.bids = writeBids(outputDir, prefix, result, R, names);
end

if any(formats == "spm")
    outputs.spm = writeSpm(outputDir, prefix, R, names);
end

if any(formats == "afni")
    outputs.afni = writeAfni(outputDir, prefix, R, names);
end

if any(formats == "conn")
    outputs.conn = writeConn(outputDir, prefix, R, names);
end

end

function file = writeMatlab(outputDir, prefix, result, R, names)
file = fullfile(outputDir, prefix + "_phase_regressors.mat");
cfg = result.config;
regressors = R;
regressorNames = cellstr(names);
TR = getOptionalField(result, 'TR', NaN);
svdDiagnostics = buildSvdDiagnostics(result);
save(file, 'regressors', 'regressorNames', 'R', 'TR', 'cfg', 'svdDiagnostics');
end

function files = writeBids(outputDir, prefix, result, R, names)
files = struct();
files.tsv = fullfile(outputDir, prefix + "_desc-fMRIPhase_timeseries.tsv");
writeTableWithHeader(files.tsv, R, names);

metadata = struct();
metadata.Description = "fMRI phase-derived nuisance regressors";
metadata.Columns = cellstr(names);
metadata.SourcePhaseFile = char(getOptionalField(result.config, 'phaseFile', ""));
metadata.SourceMagnitudeFile = char(getOptionalField(result.config, 'magnitudeFile', ""));
metadata.TR = getOptionalField(result, 'TR', NaN);
metadata.HarmonicColumns = cellstr(names(startsWith(names, "phase_harmonic_")));
metadata.SvdColumns = cellstr(names(startsWith(names, "phase_svd_")));
metadata.SvdDiagnostics = buildSvdDiagnostics(result, false);

files.json = fullfile(outputDir, prefix + "_desc-fMRIPhase_timeseries.json");
writeText(files.json, jsonencode(metadata, PrettyPrint=true));
end

function files = writeSpm(outputDir, prefix, R, names)
files = struct();

files.mat = fullfile(outputDir, prefix + "_spm_multiple_regressors.mat");
save(files.mat, 'R', 'names');

files.txt = fullfile(outputDir, prefix + "_spm_multiple_regressors.txt");
writematrix(R, files.txt, 'Delimiter', 'tab', 'FileType', 'text');
end

function file = writeAfni(outputDir, prefix, R, names)
file = fullfile(outputDir, prefix + "_afni.1D");
fid = fopen(file, 'w');
if fid < 0
    error('fmriPhase:regressors:OpenFile', 'Could not open "%s" for writing.', file);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# fMRI phase-derived nuisance regressors\n');
fprintf(fid, '# ColumnLabels = "%s"\n', strjoin(cellstr(names), ';'));
fmt = [repmat('%.10g\t', 1, size(R, 2) - 1), '%.10g\n'];
for row = 1:size(R, 1)
    fprintf(fid, fmt, R(row, :));
end
end

function files = writeConn(outputDir, prefix, R, names)
files = struct();
files.tsv = fullfile(outputDir, prefix + "_conn_covariates.tsv");
writeTableWithHeader(files.tsv, R, names);

files.mat = fullfile(outputDir, prefix + "_conn_covariates.mat");
data = R;
namesCell = cellstr(names);
save(files.mat, 'data', 'namesCell');

files.perCovariateDir = fullfile(outputDir, prefix + "_conn_covariates");
if ~isfolder(files.perCovariateDir)
    mkdir(files.perCovariateDir);
end

for idx = 1:numel(names)
    writematrix(R(:, idx), fullfile(files.perCovariateDir, names(idx) + ".txt"), ...
        'Delimiter', 'tab', 'FileType', 'text');
end
end

function writeTableWithHeader(file, R, names)
fid = fopen(file, 'w');
if fid < 0
    error('fmriPhase:regressors:OpenFile', 'Could not open "%s" for writing.', file);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strjoin(cellstr(names), sprintf('\t')));
fmt = [repmat('%.10g\t', 1, size(R, 2) - 1), '%.10g\n'];
for row = 1:size(R, 1)
    fprintf(fid, fmt, R(row, :));
end
end

function writeText(file, text)
fid = fopen(file, 'w');
if fid < 0
    error('fmriPhase:regressors:OpenFile', 'Could not open "%s" for writing.', file);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', text);
end

function value = getOptionalField(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function diagnostics = buildSvdDiagnostics(result, includeFieldMaps)
if nargin < 2
    includeFieldMaps = true;
end

diagnostics = struct();
diagnostics.SelectedComponent = NaN;
diagnostics.SelectionScores = [];
diagnostics.SingularValues = [];
diagnostics.Method = "maxDerivativeL1";
diagnostics.RespiratoryPhaseMethod = "hilbert";
diagnostics.SelectedTimecourse = [];
diagnostics.RespiratoryPhase = [];
diagnostics.RespiratoryPhaseSin = [];
diagnostics.RespiratoryPhaseCos = [];
diagnostics.LowFieldIndex = NaN;
diagnostics.HighFieldIndex = NaN;

if includeFieldMaps
    diagnostics.BreathingFieldLow = [];
    diagnostics.BreathingFieldHigh = [];
    diagnostics.BreathingFieldPeakToPeak = [];
end

if isfield(result, 'respiration')
    respiration = result.respiration;
    diagnostics.SelectedComponent = getOptionalField(respiration, 'selectedComponent', NaN);
    diagnostics.SelectionScores = getOptionalField(respiration, 'selectionScores', []);
    diagnostics.SingularValues = getOptionalField(respiration, 'singularValues', []);
    diagnostics.SelectedTimecourse = getOptionalField(respiration, 'selectedTimecourse', []);
    diagnostics.RespiratoryPhase = getOptionalField(respiration, 'respiratoryPhase', []);
    diagnostics.RespiratoryPhaseSin = getOptionalField(respiration, 'respiratoryPhaseSin', []);
    diagnostics.RespiratoryPhaseCos = getOptionalField(respiration, 'respiratoryPhaseCos', []);

    selectedTimecourse = getOptionalField(respiration, 'selectedTimecourse', []);
    respiratoryVolume = getOptionalField(respiration, 'respiratoryVolume', []);
    if ~isempty(selectedTimecourse) && ~isempty(respiratoryVolume)
        [~, lowIdx] = min(selectedTimecourse);
        [~, highIdx] = max(selectedTimecourse);
        diagnostics.LowFieldIndex = lowIdx;
        diagnostics.HighFieldIndex = highIdx;

        if includeFieldMaps
            fieldLow = single(respiratoryVolume(:, :, :, lowIdx));
            fieldHigh = single(respiratoryVolume(:, :, :, highIdx));
            diagnostics.BreathingFieldLow = fieldLow;
            diagnostics.BreathingFieldHigh = fieldHigh;
            diagnostics.BreathingFieldPeakToPeak = fieldHigh - fieldLow;
        end
    end
end
end
