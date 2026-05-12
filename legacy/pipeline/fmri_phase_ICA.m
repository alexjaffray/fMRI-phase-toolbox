%% ICA for fMRI and QSM
% 
% Independent component analysis used to see if we can separate breathing /
% cardiac

% fMRI presents the possibility of saving phase data, which can be used for
% QSM or as input to other phase-sensitive MR imaging pipelines

% 10 Aug 2022
% Editors (alphabetical order):
% Alexander Jaffray: ajaffray@physics.ubc.ca
% Michelle Medina: mmedina002@physics.ubc.ca

% PCA and ICA code adapted from: https://towardsdatascience.com/independent-component-analysis-ica-a3eba0ccec35
%

%%
clear all;
close all;

run('/home/mmedina/Desktop/fMRI_QSM/QSM/addpathqsm.m'); % change this to your own path where the QSM toolbox is stored

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
[fl, mask1] = sharp(uphas, mask, vsz,9:-2*max(vsz):2*max(vsz), 0.05);

%% Susceptibility map calculation
x = rts(fl,mask1,vsz,bdir);

%% Reshape to prepare for the SVD

p = reshape(fl,prod(imSize),s(4)); % can change between fl, x, and uphas

%% Recommended step: perform PCA on dataset since it autoscales the data
% steps taken from EEG blinking example in towardsdatascience.com

% Try different numbers
q = 21; 
% PCA function in Matlab normalizes the data (Normalization is crucial for
% ICA)
[coeff, Data_PCA, latent, tsquared, explained, mu] = pca(double(p), 'NumComponents',q);

% Explained variation used to choose q value
disp(strcat("Top ", string(q), " principle components explain ", string(sum(explained(1:q))), " of variation"));

%% Use ICA
% Matlab only has the function RICA; based on ICA byt has a reconstruction
% cost (penalty or regulatization term)

% Principal components will be used to compute the independent components
Mdl = rica(Data_PCA, q);

% ICA application
Data_ICA = transform(Mdl, Data_PCA);

%% Plot Components

% Number of plots/column 
plts = 7;

% Plot
figure(2)
for i=1:q
    subplot(plts, ceil(q/plts),i)
    plot(Data_ICA(:,i).^2)
    title(strcat("Component ", string(i), " Squared"))
    ax = gca;
    ax.XTickLabel = {};
end


