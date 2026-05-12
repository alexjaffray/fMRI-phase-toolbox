function startup_fmri_phase()
%STARTUP_FMRI_PHASE Add the v2 MATLAB toolbox to the MATLAB path.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir, 'matlab'));
addpath(fullfile(rootDir, 'matlab', 'vendor', 'tapas_physio'));

end
