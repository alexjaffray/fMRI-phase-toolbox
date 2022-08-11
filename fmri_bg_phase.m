%% Pipeline for processing fMRI phase and magnitude data for QSM using the QSM toolbox
% 
% fMRI presents the possibility of saving phase data, which can be used for
% QSM or as input to other phase-sensitive MR imaging pipelines

% 8 Aug 2022
% Editors (alphabetical order):
% Alexander Jaffray: ajaffray@physics.ubc.ca
%
%

%%
clear all;
close all;

run('/srv/data/ajaffray/QSM/addpathqsm.m'); % change this to your own path where the QSM toolbox is stored

%% Get magnitude and phase data from the nifti files
angleData = niftiread(uigetfile());
magnitudeData = niftiread(uigetfile());

%% Set these based on the parameters of the imaging experiment
TE = 0.03; %[s]
B0 = 3.0; %[T]
GYRO = 267.513; %[rad/s/uT]

%% Voxel size, b0 direction and image size
vsz = [2.4962 2.4962 2.5000]; %[mm]
bdir = [0 0 1]; %[arbitrary]
s = size(angleData);
imSize = [s(1),s(2),s(3)]; %[voxels]

%% Generate the mask from the last timepoint magnitude
mask = generateMask(magnitudeData(:,:,:,end), vsz, '-m -n -f 0.5');

%% Convert from Int16 to phase 
phas = double(angleData) ./ 2048 * pi - pi;

%% Do the laplacian unwrapping and rescale to units of field
uphas = unwrapLaplacian(phas,mask,vsz);
uphas = uphas ./ (B0 * GYRO * TE);

%% Background field removal
[fl, mask1] = resharp(uphas, mask, vsz,9:-2*max(vsz):2*max(vsz), 0.05);

%% Susceptibility map calculation
x = rts(fl,mask1,vsz,bdir);

%% Reshape to prepare for the SVD
p = reshape(fl,prod(imSize),s(4));

%% Decompose the X Map using the SVD
[U,S,V] = svd(double(p),"econ");

%% Take the first 5 components of the SVD

% first component
componentVector = 1;
uphas1 = recomposeSVD(U,S,V,componentVector,imSize);

% second component
componentVector = 2;
uphas2 = recomposeSVD(U,S,V,componentVector,imSize);

% third component
componentVector = 3;
uphas3 = recomposeSVD(U,S,V,componentVector,imSize);


% fourth component
componentVector = 4;
uphas4 = recomposeSVD(U,S,V,componentVector,imSize);


% fifth component
componentVector = 5;
uphas5 = recomposeSVD(U,S,V,componentVector,imSize);


%% Interpolate whole volume using zero-filling along the time dimension to help visualize things

interpolationFactor = 3;

% interpolate the original X map
interpAbs = interpft(fl,260*interpolationFactor,4); 

% interpolate the 5 SVD components
interpTime1 = interpft(uphas1,260*interpolationFactor,4);
interpTime2 = interpft(uphas2,260*interpolationFactor,4);
interpTime3 = interpft(uphas3,260*interpolationFactor,4);
interpTime4 = interpft(uphas4,260*interpolationFactor,4);
interpTime5 = interpft(uphas5,260*interpolationFactor,4);

%% Create the time vector from the interpolated data

timeVector = 0:1.15/interpolationFactor:260*1.15
timeVector = timeVector(1:260*interpolationFactor);


%% Select a Slice
selectedSlice = 23;

%% Show the variation of the data in real time by dynamic plotting or single plotting

doDynamicPlot = true;

if doDynamicPlot

    figure(22);
    set(gcf,'Position',[100 100 1500 1000])
    %vidfile = VideoWriter('testmovieCHI.avi');
    %open(vidfile);

    for ind = 1:260*interpolationFactor
        subplot(2,3,1);
        im = squeeze(interpAbs(:,:,selectedSlice,ind));
        imagesc(im),colormap(gray),caxis([-0.2,0.2]),colorbar;
        title('Original X Image');
        subplot(2,3,2);
        im2 = squeeze(interpTime1(:,:,selectedSlice,ind));
        imagesc(im2),colormap(gray),caxis([-0.2,0.2]),colorbar;
        title('SVD 1');
        subplot(2,3,3);
        im3 = squeeze(interpTime2(:,:,selectedSlice,ind));
        imagesc(im3),colormap(gray),caxis([-0.15,0.15]),colorbar; 
        title('SVD 2');
        subplot(2,3,4);
        im4 = squeeze(interpTime3(:,:,selectedSlice,ind));
        imagesc(im4),colormap(gray),caxis([-0.1,0.1]),colorbar;
        title('SVD 3');
        subplot(2,3,5);
        im5 = squeeze(interpTime4(:,:,selectedSlice,ind));
        imagesc(im5),colormap(gray),caxis([-0.1,0.1]),colorbar;
        title('SVD 4');
        subplot(2,3,6);
        im6 = squeeze(interpTime5(:,:,selectedSlice,ind));
        imagesc(im6),colormap(gray),caxis([-0.1,0.1]),colorbar;
        title('SVD 5');
        drawnow

        %F(ind) = getframe(gcf);
        %writeVideo(vidfile,F(ind));
    end
    %close(vidfile);

end

%% Define Ranges in the slice to look at Fluctuation
range1 = 46:48;
range2 = 26:28;

range1 = 23:31;
range2 = 25:45;

%% Plot in real time the first 5 SVD coefficients in the region of interest defined by range2 and range1

doDynamicPlot = false;

if doDynamicPlot

    figure(19);
    set(gcf,'Position',[100 100 1500 1000])
    title('Field Fluctuations: Original Image and First 5 SVD Coefficients');
    xlabel('Time (s)');
    ylabel('Field Fluctuations (ppm)');
    %vidfile = VideoWriter('testmovieCHI_FLUX.avi');
    %open(vidfile);

    h1 = animatedline('Color','black');
    h2 = animatedline('Color','red');
    h3 = animatedline('Color','blue');
    h4 = animatedline('Color','black','LineStyle','-.');
    h5 = animatedline('Color','cyan');
    h6 = animatedline('Color','magenta');
    axis([0 1.15*259 -0.005 0.005]);

    for ind = 1:260*3
        addpoints(h1,timeVector(ind),mean(interpTime1(range2,range1,selectedSlice,ind),'all'));
        drawnow
        addpoints(h2,timeVector(ind),mean(interpTime2(range2,range1,selectedSlice,ind),'all'));
        drawnow
        addpoints(h3,timeVector(ind),mean(interpTime3(range2,range1,selectedSlice,ind),'all'));
        drawnow
        addpoints(h5,timeVector(ind),mean(interpTime4(range2,range1,selectedSlice,ind),'all'));
        drawnow
        addpoints(h6,timeVector(ind),mean(interpTime5(range2,range1,selectedSlice,ind),'all'));
        drawnow
        addpoints(h4,timeVector(ind),mean(interpAbs(range2,range1,selectedSlice,ind),'all'));
        drawnow
        legend('1st Component','2nd Component','3rd Component','4th Component','5th Component','Original Image');
        %F(ind) = getframe(gcf);
        %writeVideo(vidfile,F(ind));
    end
    %close(vidfile);

else
    
    figure(20);
    set(gcf,'Position',[100 100 1500 1000])
    title('Field Fluctuations: Original Image and First 5 SVD Coefficients');
    xlabel('Time (s)');
    ylabel('Field Fluctuations (ppm)');
    %vidfile = VideoWriter('testmovieCHI_FLUX.avi');
    %open(vidfile);
    
    plot(timeVector,squeeze(mean(interpAbs(range2,range1,selectedSlice,:),[1 2 3])),'Color','black','LineStyle','-.');
    hold on;
    plot(timeVector,squeeze(mean(interpTime1(range2,range1,selectedSlice,:),[1 2 3])),'Color','red');
    hold on;
    plot(timeVector,squeeze(mean(interpTime2(range2,range1,selectedSlice,:),[1 2 3])),'Color','blue');
    hold on;
    plot(timeVector,squeeze(mean(interpTime3(range2,range1,selectedSlice,:),[1 2 3])),'Color','black');
    hold on;
    plot(timeVector,squeeze(mean(interpTime4(range2,range1,selectedSlice,:),[1 2 3])),'Color','cyan');
    hold on;
    plot(timeVector,squeeze(mean(interpTime5(range2,range1,selectedSlice,:),[1 2 3])),'Color','magenta');

    axis([0 1.15*259 -0.5 0.5]);
    
    legend('Original Image','1st Component','2nd Component','3rd Component','4th Component','5th Component');

end

figure();
set(gcf,'Position',[100 100 1500 1000])
title('SVD coefficients');
xlabel('Singular Value Number');
ylabel('Singular Value Coefficient');
plot(1:s(4),S*ones(s(4),1));
set(gca,'yscale','log');
    

%% Extract out the 2nd SVD component from the region of interest defined by range2 and range1 and take it's spatial mean
d = squeeze(mean(interpTime2(range2,range1,selectedSlice,:),[1 2 3]));

%% Plot the time-course of the 2nd SVD coefficient and identify minima

% identify the minima (for gating of breath for example)
foundMins = islocalmin(d,'MinProminence',0.0035);

% Plot
figure();
plot(timeVector,d);
hold on;
scatter(timeVector(foundMins),d(foundMins),'o');
title('2nd SVD Coefficient');
xlabel('Time (s)');
ylabel('Fluctuation in X map (ppm)');
legend('Delta X','Minima');

%% Recompose the image from the Nth component of the SVD
function recomposed = recomposeSVD(U,S,V,coeffVec,msize)
    mask1 = zeros(size(S));
    mask1(coeffVec,coeffVec) = 1;
    S = S .* mask1;
    recomposed = reshape(U*S*V',[msize,size(U,2)]);

end