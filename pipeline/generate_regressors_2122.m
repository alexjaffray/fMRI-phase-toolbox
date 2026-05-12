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
%%
baseDir = '/data/ubcitar/ICE_analysis/raw/';
outputDir = '/srv/data/ajaffray/concussion_study/regressors';

%% 2021-2022
subjects = {'sub-M01','sub-M02','sub-M03','sub-M04','sub-M05','sub-M06','sub-M07','sub-M08','sub-M09','sub-M10','sub-M11','sub-M13','sub-M14','sub-M15','sub-M16','sub-M17','sub-M18','sub-M19','sub-M20','sub-M21','sub-M23','sub-M24','sub-M25','sub-M26','sub-M29','sub-MC02','sub-MC03','sub-MC05','sub-MC06','sub-MC08','sub-MC09','sub-MC10','sub-MC11','sub-MC14','sub-MC15','sub-W01','sub-W02','sub-W03','sub-W04','sub-W05','sub-W06','sub-W07','sub-WC01','sub-WC02','sub-WC03','sub-WC04','sub-WC05','sub-WC06','sub-WC08','sub-WC09','sub-WC11'};
seasons = {'ses-2122begin','ses-2122end'};

nsub = length(subjects);
nseas = length(seasons);

for ii = 1:length(subjects)
    for jj = 1:length(seasons)
        subject_list{(jj-1)*nsub + ii} = strjoin({subjects{ii},seasons{jj}},'_');
        subject_path{(jj-1)*nsub + ii} = fullfile(subjects{ii},seasons{jj});
    end
end

%%
for ii = 1:length(subject_list)
    % Set these for current protocol

    subprefix = subject_list{ii}

    if ~exist(fullfile(outputDir,strjoin({subprefix,'reg.mat'},'_')))
    
    angleFile = strjoin({subprefix,'task-rest_part-phase_bold.nii'},'_');
    magFile = strjoin({subprefix,'task-rest_part-mag_bold.nii'},'_');

    datDir = fullfile(baseDir,subject_path{ii},'func');
    
    % if ~exist(angleFile)
    % 
    %     %% Read in scan info from template NIFTI
    %     [angleFile,angleDir] = uigetfile("*.nii","Select the Template Phase NIFTI File");
    %     [magFile,magDir] = uigetfile("*.nii","Select the Template Magnitude NIFTI File");   
    % end

    if exist(fullfile(datDir,angleFile))

        % Get magnitude and phase data from the nifti files
        angleData = niftiread(fullfile(datDir,angleFile));
        magnitudeData = niftiread(fullfile(datDir,magFile));
    
        % doRotate
        angleData = rot90(angleData);
        magnitudeData = rot90(magnitudeData);
        
        % Set scan info
        scaninfo = niftiinfo(fullfile(datDir,angleFile));
        
        % Set these based on the parameters of the imaging experiment (HARDCODED for now!)
        TE = 0.03; %[s]
        B0 = 3.0; %[T]
        GYRO = 267.513; %[rad/s/uT]
        
        % Voxel size, b0 direction and image size
        vsz = scaninfo.PixelDimensions(1:3); %[mm]
        bdir = [0 0 1]; %[arbitrary]
        s = scaninfo.ImageSize;
        imSize = [s(1),s(2),s(3)]; %[voxels]
        TR = scaninfo.PixelDimensions(4);
        
        % Generate the mask from the last timepoint magnitude
        mask0 = generateMask(magnitudeData(:,:,:,end), vsz, '-m -n -f 0.5');
        
        % Generate harmonic fields
        [harmfields,mask1,phas,uphas,fl] = calcHarmFields(angleData,mask0,TE,B0,GYRO,vsz);
        
        % Choose inputImage and Reshape to prepare for the SVD
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
        
        % Set SVD params
        interpolationFactor = 1;
        exampleSlice = 31;
        doPlot = true;
        
        plotSlice(phas,uphas,fl,harmfields,exampleSlice);
        
        % Decompose the Image using the SVD
        [respvol,timeVector,zero_ord_vol, other_comps] = getRespComp(inputImage,s,TR,interpolationFactor,doPlot);
        
        % Test Slice TR
        % [respvol,timeVector,rawresp] = getRespCompSliceTR(inputImage,s,TR,interpolationFactor);
        
        % Project Respiratory correlated field onto spherical harmonics
        
        % define a regular grid covering the domain of the data (arb scale)
        [xxmat,yymat,zzmat] = meshgrid(-47.5:1:47.5,-47.5:1:47.5,-28:1:28);
        order = ones(96,96,57);
        
        circmask = ones(96,96,57);
        circmask = circmask((xxmat.^2 + yymat.^2 + zzmat.^2).^0.5 < 40);
        
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
        
        % Ignore slice timing -> seems more promising :) 
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
        
        % Create Regressor Matrix
        regmat.reg = normalize(cat(2, volTR_coeffs, other_comps),1);
        save(fullfile(outputDir,strjoin({subprefix,'reg.mat'},'_')),'-struct','regmat')

    else
        
        subprefix

    end

    end

    ii

end

