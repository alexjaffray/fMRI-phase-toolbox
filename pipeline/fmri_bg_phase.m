%% Pipeline for processing fMRI phase and magnitude data for QSM using the QSM toolbox
% 
% fMRI presents the possibility of saving phase data, which can be used for
% QSM or as input to other phase-sensitive MR imaging pipelines

% 8 Aug 2022
% Editors (alphabetical order):
% Alexander Jaffray: ajaffray@physics.ubc.ca
% Michelle Medina: mmedina002@phas.ubc.ca
%

%%
clear all;
close all;

%%

run('/srv/data/ajaffray/QSM/addpathqsm.m'); % change this to your own path where the QSM toolbox is stored

%%
run('/srv/data/ajaffray/MRecon-5.0.11/startup.m');

%% Get magnitude and phase data from the nifti files

[angleFile,angleDir] = uigetfile("*.nii","Select the Phase NIFTI File");
[magFile,magDir] = uigetfile("*.nii","Select the Magnitude NIFTI File");

%%
angleData = niftiread(fullfile(angleDir,angleFile));
magnitudeData = niftiread(fullfile(magDir,magFile));

%% data From MRecon
mreconDat = MRecon();
mreconDat.ReadData();


%%
angleData = squeeze(mreconDat.Data(:,:,:,1,:,1,2));
magnitudeData = squeeze(mreconDat.Data(:,:,:,1,:,1,1));

%% Read in scan info

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

mask = generateMask(magnitudeData(:,:,:,end), vsz, '-m -n -f 0.5');

%% Convert from Int16 to phase 

phas = double(angleData) ./ 2048 * pi - pi;

%% Do the laplacian unwrapping and rescale to units of field

uphas = unwrapLaplacian(phas,mask,vsz);
uphas = uphas ./ (B0 * GYRO * TE);

%% Background field removal

[fl, mask1] = resharp(uphas, mask, vsz,9:-2*max(vsz):2*max(vsz), 0.05);

%% get the harmonic phase evolution

harmfields = uphas - fl;

%% Susceptibility map calculation

x = rts(fl,mask1,vsz,bdir);

%% Choose inputImage and Reshape to prepare for the SVD
inputImage = zeros(size(uphas));

imageType = "totalField";

switch imageType
    case "chi"
        inputImage = x;
    case "totalField"
        inputImage = uphas;
    case "localField"
        inputImage = fl;
end

p = reshape(inputImage,prod(imSize),s(4));

%% Decompose the Image using the SVD

[U,S,V] = svd(double(p),"econ");

% Take the first 5 components of the SVD

% first component
componentVector = 1;
comp1 = recomposeSVD(U,S,V,componentVector,imSize);

% second component
componentVector = 2;
comp2 = recomposeSVD(U,S,V,componentVector,imSize);

% third component
componentVector = 3;
comp3 = recomposeSVD(U,S,V,componentVector,imSize);

% fourth component
componentVector = 4;
comp4 = recomposeSVD(U,S,V,componentVector,imSize);

% fifth component
componentVector = 5;
comp5 = recomposeSVD(U,S,V,componentVector,imSize);

%% Interpolate whole volume using zero-filling along the time dimension to help visualize things
doInterpolation = true;

interpAbs = [];
interpTime1 = [];
interpTime2 = [];
interpTime3 = [];
interpTime4 = [];
interpTime5 = [];

interpolationFactor = 1;

if doInterpolation

    interpolationFactor = 2;

    % interpolate the original image in time
    interpAbs = interpft(inputImage,s(4)*interpolationFactor,4); 

    % interpolate the 5 SVD components in time
    interpTime1 = interpft(comp1,s(4)*interpolationFactor,4);
    interpTime2 = interpft(comp2,s(4)*interpolationFactor,4);
    interpTime3 = interpft(comp3,s(4)*interpolationFactor,4);
    interpTime4 = interpft(comp4,s(4)*interpolationFactor,4);
    interpTime5 = interpft(comp5,s(4)*interpolationFactor,4);

else
   
    interpAbs = inputImage;
    interpTime1 = comp1;
    interpTime2 = comp2;
    interpTime3 = comp3;
    interpTime4 = comp4;
    interpTime5 = comp5;
end

%% Create the time vector from the interpolated data

timeVector = 0:TR/interpolationFactor:s(4)*TR;
timeVector = timeVector(1:s(4)*interpolationFactor);


%% Define Ranges in the slice to look at Fluctuation and calculate roi means

selectedSlice = 31;

width = 8;
height = 20;
xmin = 40;
ymin = 40;

range2 = xmin:(xmin+height);
range1 = ymin:(ymin+width);


d0 = getROImean(interpAbs,range2,range1,selectedSlice,TR/interpolationFactor);
d1 = getROImean(interpTime1,range2,range1,selectedSlice,TR/interpolationFactor);
d2 = getROImean(interpTime2,range2,range1,selectedSlice,TR/interpolationFactor);
d3 = getROImean(interpTime3,range2,range1,selectedSlice,TR/interpolationFactor);
d4 = getROImean(interpTime4,range2,range1,selectedSlice,TR/interpolationFactor);
d5 = getROImean(interpTime5,range2,range1,selectedSlice,TR/interpolationFactor);

% %% test ICA!
% 
% % define data
% 
% XDat = inputImage(range2,range1,selectedSlice,:);
% XDatReshaped = reshape(XDat,[],size(XDat,4));
% 
% model = rica(XDatReshaped',6);
% 
% reconstructedData = transform(model, XDatReshaped');
% 
% figure();
% plot(reconstructedData);

%% Define axesranges for the plotting


dataMin = min(d0,[],'all');
dataMax = max(d0,[],'all');
dataRange = dataMax - dataMin;
dlim = [0 0];

if dataMin < 0
    
    dlim(1) = dataMin - dataRange;
    dlim(2) = 2*dataRange;
    
else
    dlim(1) = -2*dataRange;
    dlim(2) = dataMax + dataRange;
end


%% Plot the first 5 SVD coefficients in the region of interest defined by range2 and range1

close all;

figure();
set(gcf,'Position',[100 100 2000 750])
subplot(1,2,1);
imagesc(interpAbs(:,:,selectedSlice,1));
title('Image with ROI');
xlabel('x-axis');
ylabel('y-axis');
colormap(gray);
ax = gca;
hold on;

roi = drawrectangle(ax,'Position',[xmin,ymin,width,height]);

subplot(1,2,2);
plot(timeVector,d0,'Color','black','LineStyle','-.');
hold on;
plot(timeVector,d1,'Color','red');
hold on;
plot(timeVector,d2,'Color','blue');
hold on;
plot(timeVector,d3,'Color','black');
hold on;
plot(timeVector,d4,'Color','cyan');
hold on;
plot(timeVector,d5,'Color','magenta');

title('Field Fluctuations: Original Image and First 5 SVD Coefficients');
xlabel('Time (s)');
ylabel('Field Fluctuations (ppm)');

axis([timeVector(1) timeVector(end) dlim(1) dlim(2)]);

legend('Original Image','1st Component','2nd Component','3rd Component','4th Component','5th Component');

%% Show the variation of the data in real time by dynamic plotting or single plotting

doDynamicPlot = false;

if doDynamicPlot

    figure(22);
    set(gcf,'Position',[100 100 1500 1000])
    %vidfile = VideoWriter('testmovieCHI.avi');
    %open(vidfile);

    for ind = 1:s(4)*interpolationFactor
        subplot(2,3,1);
        im = squeeze(interpAbs(:,:,selectedSlice,ind));
        imagesc(im),colormap(gray),caxis([-0.05,0.05]),colorbar;
        title('Original X Image');
        subplot(2,3,2);
        im2 = squeeze(interpTime1(:,:,selectedSlice,ind));
        imagesc(im2),colormap(gray),caxis([-0.05,0.05]),colorbar;
        title('SVD 1');
        subplot(2,3,3);
        im3 = squeeze(interpTime2(:,:,selectedSlice,ind));
        imagesc(im3),colormap(gray),caxis([-0.015,0.015]),colorbar; 
        title('SVD 2');
        subplot(2,3,4);
        im4 = squeeze(interpTime3(:,:,selectedSlice,ind));
        imagesc(im4),colormap(gray),caxis([-0.01,0.01]),colorbar;
        title('SVD 3');
        subplot(2,3,5);
        im5 = squeeze(interpTime4(:,:,selectedSlice,ind));
        imagesc(im5),colormap(gray),caxis([-0.01,0.01]),colorbar;
        title('SVD 4');
        subplot(2,3,6);
        im6 = squeeze(interpTime5(:,:,selectedSlice,ind));
        imagesc(im6),colormap(gray),caxis([-0.01,0.01]),colorbar;
        title('SVD 5');
        drawnow

        %F(ind) = getframe(gcf);
        %writeVideo(vidfile,F(ind));
    end
    %close(vidfile);

end


%% Plot the singular values

figure();
set(gcf,'Position',[100 100 1500 1000]);
plot(1:s(4),S*ones(s(4),1));
title('Singular Values');
xlabel('Singular Value Number');
ylabel('Singular Value');
set(gca,'yscale','log');

%% Plot the time-course of the 2nd SVD coefficient and identify minima

% identify the minima (for gating of breath for example)
foundMins = islocalmin(d3,'MinProminence',0.0035);

% Plot
figure();
plot(timeVector,d3);
hold on;
scatter(timeVector(foundMins),d3(foundMins),'o');
title('2nd SVD Coefficient');
xlabel('Time (s)');
ylabel('Fluctuation in X map (ppm)');
legend('Delta X','Minima');

%% Plot the standard deviation of the qsm values obtained throughout the acquisition

xRange = std(uphas,0,4);
figure();
imagesc(xRange(:,:,25)),colormap('hot'),caxis([0 0.05]);
colorbar;
