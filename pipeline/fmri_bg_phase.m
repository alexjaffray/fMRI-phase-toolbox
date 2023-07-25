%% Pipeline for processing fMRI phase and magnitude data for QSM using the QSM toolbox
%
% fMRI presents the possibility of saving phase data, which can be used for
% QSM or as input to other phase-sensitive MR imaging pipelines

% 29 Nov 2022
% Editors (alphabetical order):
% Alexander Jaffray: ajaffray@physics.ubc.ca
% Michelle Medina: mmedina002@phas.ubc.ca  
clear all;
close all;

%%
addpath('/srv/data/ajaffray/fMRI-phase-toolbox/');

%%
run('/srv/data/ajaffray/QSM/addpathqsm.m'); % change this to your own path where the QSM toolbox is stored

%%
run('/srv/data/ajaffray/MRecon-5.0.11/startup.m');

%% Change this for each run
fileLabel = "/srv/data/ajaffray/fMRI-phase-toolbox/data/christina.mat";

%% Set these for current protocol
angleFile = 'sub-M02_ses-2122post_task-rest_run-1_part-phase_bold.nii';
angleDir = '/srv/data/scratch/sub-M02/ses-2122post/func/';
magDir = '/srv/data/scratch/sub-M02/ses-2122post/func/';
magFile = 'sub-M02_ses-2122post_task-rest_run-1_part-mag_bold.nii';
% 
% if ~exist(angleFile)
% 
%     %% Read in scan info from template NIFTI
%     [angleFile,angleDir] = uigetfile("*.nii","Select the Template Phase NIFTI File");
%     [magFile,magDir] = uigetfile("*.nii","Select the Template Magnitude NIFTI File");   
% end

%% Ask for data files (.rec, .par files must have lowercase file tails)
dataFormat = "rec";

switch dataFormat
    case "rec"
        % data From MRecon
        % filename has to be in lower case!!!
        mreconDat = MRecon();
        mreconDat.ReadData();
        angleData = squeeze(mreconDat.Data(:,:,:,1,:,1,2));
        magnitudeData = squeeze(mreconDat.Data(:,:,:,1,:,1,1));
        
    case "nifti"
        % Get magnitude and phase data from the nifti files
        [angleFile,angleDir] = uigetfile("*.nii","Select the Phase NIFTI File");
        [magFile,magDir] = uigetfile("*.nii","Select the Magnitude NIFTI File");
        angleData = niftiread(fullfile(angleDir,angleFile));
        magnitudeData = niftiread(fullfile(magDir,magFile));

end

%% Ask User for the phys log file!
[logFile,logDir] = uigetfile("*.log","Select the Relevant PhysLog File");
physLogTable = readPhysLog(fullfile(logDir,logFile));

%% Set scan info
scaninfo = niftiinfo(fullfile(angleDir,angleFile));

%% Set these based on the parameters of the imaging experiment (HARDCODED for now!)
TE = 0.03; %[s]
B0 = 3.0; %[T]
GYRO = 267.513; %[rad/s/uT]

%% Voxel size, b0 direction and image size
vsz = scaninfo.PixelDimensions(1:3); %[mm]
bdir = [0 0 1]; %[arbitrary]
s = scaninfo.ImageSize;
imSize = [s(1),s(2),s(3)]; %[voxels]
TR = scaninfo.PixelDimensions(4);

%% Generate the mask from the last timepoint magnitude
mask0 = generateMask(magnitudeData(:,:,:,end), vsz, '-m -n -f 0.5');

%% Generate harmonic fields
[harmfields,mask1,phas,uphas,fl] = calcHarmFields(angleData,mask0,TE,B0,GYRO,vsz);

%% Choose inputImage and Reshape to prepare for the SVD
inputImage = zeros(size(uphas));

imageType = "harmonicField";

switch imageType
    case "chi"
        % Susceptibility map calculation
        x = rts(fl,mask1,vsz,bdir);
        inputImage = x;
    case "totalField"
        inputImage = uphas;
    case "localField"
        inputImage = fl;
    case "harmonicField"
        inputImage = harmfields;
end

%% Set SVD params
interpolationFactor = 1;
exampleSlice = 31;
doPlot = true;

plotSlice(phas,uphas,fl,harmfields,exampleSlice);

%% Decompose the Image using the SVD
[respvol,timeVector,zero_ord_vol] = getRespComp(inputImage,s,TR,interpolationFactor,doPlot);

%% Project Respiratory correlated field onto spherical harmonics

% define a regular grid covering the domain of the data (arb scale)
[xxmat,yymat,zzmat] = meshgrid(-47.5:1:47.5,-47.5:1:47.5,-28:1:28);
order = ones(96,96,57);

circmask = ones(96,96,57);
circmask = circmask((xxmat.^2 + yymat.^2 + zzmat.^2).^0.5 .< 40)

n_timepoints = length(timeVector);
mb = 3;
n_exc = size(order,3)/mb;

xx0 = xxmat(:,:,1);
xx = xx0(:);
yy0 = yymat(:,:,1);
yy = yy0(:);

slicetrack = repmat(1:19,1,n_timepoints)';
se = strel('cube',7);

mask = imerode(mask0,se);

C = zeros(n_exc,n_timepoints,16);

elvec = 1:96;

dmat = [];

tic
for volidx = 1:n_timepoints
    for ind = 1:n_exc

        zvec = [ind ind+19 ind+38];

        xx = xxmat(mask(elvec,elvec,zvec));
        yy = yymat(mask(elvec,elvec,zvec));
        zz = zzmat(mask(elvec,elvec,zvec));

        onevec = order(mask(elvec,elvec,zvec));

        dmat = [onevec xx yy zz xx.*yy yy.*zz xx.*zz xx.^2 - yy.^2 2*zz.^2 - xx.^2 - yy.^2 xx.*yy.*zz zz.*xx.^2 - zz.*yy.^2 3*yy.*xx.^2 - yy.^3 (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*yy (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*xx 5*zz.^3-3*(xx.^2 + yy.^2 + zz.^2).*zz xx.^3 - 3*xx.*yy.^2];

        C(ind,volidx,:) = basisExpansion(respvol(elvec,elvec,zvec,volidx),dmat,mask(elvec,elvec,zvec));
    end
end
toc

C3 = reshape(C,[],16);

% Filtering
fs = n_exc/TR;
f0 = 1/TR;
fn = fs/2;
freqRatio = f0/fn;

notchWidth = 0.1;

% Compute zeros
notchZeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

% Compute poles
notchPoles = (1-notchWidth) * notchZeros;

b_notch = poly( notchZeros ); %  Get moving average filter coefficients
a_notch = poly( notchPoles ); %  Get autoregressive filter coefficients

% filter signal x
C4 = filter(b_notch,a_notch,C3);
[b_lp,a_lp] = butter(12,freqRatio/2);
sliceTR_coeffs = filter(b_lp,a_lp,C4);

%% Ignore slice timing -> seems more promising :) 
testharmfields2 = respvol;
xx = xxmat(mask);
yy = yymat(mask);
zz = zzmat(mask);

onevec = order(mask);

dmat = [onevec xx yy zz xx.*yy yy.*zz xx.*zz xx.^2 - yy.^2 2*zz.^2 - xx.^2 - yy.^2 xx.*yy.*zz zz.*xx.^2 - zz.*yy.^2 3*yy.*xx.^2 - yy.^3 (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*yy (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*xx 5*zz.^3-3*(xx.^2 + yy.^2 + zz.^2).*zz xx.^3 - 3*xx.*yy.^2];

volTR_coeffs = zeros(n_timepoints,16);

tic
for volidx = 1:n_timepoints

    volTR_coeffs(volidx,:) = basisExpansion(respvol(:,:,:,volidx),dmat,mask);

end

toc

%% Plot Each Coefficient of Varying Order (volume TR)
plotSphericalHarmonics(volTR_coeffs,TR/interpolationFactor*(1:n_timepoints));

%% Plot Each Coefficient of Varying Order (Slice TR, Filtered Notch + LP)
plotSphericalHarmonics(sliceTR_coeffs,TR/interpolationFactor/n_exc*(1:(n_timepoints*n_exc)));

%% respcomp choice
respcompVOL = volTR_coeffs(:,1); % NEED
respcompSL = sliceTR_coeffs(:,1); % NEED

%% Prepare physLog
pLogFileName = string(fullfile(logDir,logFile));
logfiles.cardiac = pLogFileName;
logfiles.respiration = pLogFileName;
logfiles.sampling_interval = 1/500;
logfiles.relative_start_acquisition = 0;
phaseSign = 1; % NEED
doTRComp = 0; % NEED

if phaseSign == -1
    doTRComp = 1;
end

[c, r, t2, cpulse, acq_codes] = tapas_physio_read_physlogfiles_philips(logfiles, 'PPU');

c = c(1:2:end);
r = r(1:2:end);
t2 = t2(1:2:end)/2;
scanStart = find(physLogTable.mark>20);
t2 = t2(scanStart+1:end);
[physLogResp,fh] = tapas_physio_filter_respiratory(r(scanStart+1:end),t2(2)-t2(1),[],true,true); 



%% SLICE TR Processing
TR = 1.150/19;

timeVector = (1:length(respcompSL)) * 1.15/19*interpolationFactor - 1.15/19*interpolationFactor;
scanTime = mreconDat.Parameter.Labels.ScanDuration;
phaseTimeSL = timeVector(end);
samplingRate = length(physLogResp) / scanTime;

sliceTime = timeVector + (scanTime - phaseTimeSL + 0*TR/2)/2 + doTRComp*TR/2 * 19; % NEED
physLogTime = linspace(0,scanTime,length(physLogResp)); % NEED

timeVectorVol = (1:length(respcompVOL))*1.15*interpolationFactor - 1.15*interpolationFactor;
phaseTimeVol = timeVectorVol(end);
volTime = timeVectorVol + (scanTime - phaseTimeVol) + doTRComp*1.15/19; % NEED

%%
% % Plot the processed respiratory phase (naive approach)
% filteredTrace = lowpass(respcompSL,0.1,interpolationFactor/TR); % low pass filter to get the jumps out of the data
% respPhase = calculateRespPhase(respcompSL);
% figure();
% plot(fmriTime,respPhase);
% title('Respiratory Phase during FMRI Acquisition');
% xlabel('Time (s)');
% ylabel('Respiratory Phase \in [-\pi, \pi]');
% hold on;
% % Load the physlog file and process as before (naive approach)
% resampledPhysLogTrace = physLogResp; %resample(physLogResp,interpolationFactor,round(samplingRate*TR)); % resample to same rate as the fmri-derived resp data
% filteredPhysLogTrace = lowpass(resampledPhysLogTrace,20,520); % low pass
% filter to remove weird data jumps
% respPhasePhysLog = calculateRespPhase(filteredPhysLogTrace);
% plot(physLogTime,respPhasePhysLog);
% xlabel('Time (s)');
% ylabel('Respiratory Phase \in [-\pi, \pi]');
% legend('fMRI phase data (Slice TR)','breathing belt data');

%% Process Vol TR Phase Data



%% Generate Figure
% The standard values for colors saved in PLOT_STANDARDS() will be accessed from the variable PS
close all;

PS = PLOT_STANDARDS();

figure(1);
fig1_comps.fig = gcf;
hold on
fig1_comps.p1 = plot(volTime,respcompVOL./max(respcompVOL));
fig1_comps.p2 = plot(sliceTime,phaseSign*respcompSL./max(respcompSL));
fig1_comps.p3 = plot(physLogTime,physLogResp);

% ADD LABELS, TITLE, LEGEND
title('Respiratory Trace during FMRI Acquisition');
xlabel('Time (s)');
ylabel('Respiratory Trace (a.u)');
legend([fig1_comps.p1, fig1_comps.p2, fig1_comps.p3], 'B_0 Fluctuation (Volume TR)','B_0 Fluctuation (Slice TR)','Breathing Belt');
legendX = .82; legendY = .87; legendWidth = 0.02; legendHeight = 0.02;
fig1_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];
% If you want the tightest box set width and height values very low matlab automatically sets the tightest box
%========================================================
% SET PLOT PROPERTIES
% Choices for COLORS can be found in ColorPalette.png
set(fig1_comps.p1, 'LineStyle', '-', 'LineWidth', 2, 'Color',PS.Green1);
set(fig1_comps.p2, 'LineStyle', '-', 'LineWidth', 2, 'Color',PS.Blue3);
set(fig1_comps.p3, 'LineStyle', '-.', 'LineWidth', 2, 'Color', PS.MyRed);
%========================================================
% INSTANTLY IMPROVE AESTHETICS-most important step
STANDARDIZE_FIGURE(fig1_comps);

% %% Calculate the peak to peak variation between two resp traces
% ts1 = timeseries(respPhase,fmriTime);
% ts2 = timeseries(respPhasePhysLog,physLogTime);
% 
% [tsout1,tsout2] = synchronize(ts1,ts2,'union');
% 
% [peaks,locs,w] = findpeaks(tsout1.Data);
% [peaks2,locs2,w2] = findpeaks(tsout2.Data);
% 
% idx2 = [1:22 24:49];
% 
% rms(tsout1.Time(locs(idx2)) - tsout2.Time(locs2))

%% Save Things for Plotting
save(fileLabel,"sliceTR_coeffs","volTR_coeffs","sliceTime","volTime","physLogTime","physLogResp","respcompVOL","respcompSL","phaseSign", "doTRComp");

%% Save the mask
save("/srv/data/ajaffray/fMRI-phase-toolbox/data/mask_christina.mat","mask0");