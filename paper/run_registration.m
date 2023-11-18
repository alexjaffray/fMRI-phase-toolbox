close all;
clear all;

load imdatamat.mat

%%
fixed = magnitudeData(:,:,:,1);
moving = magnitudeData(:,:,:,2);

%% View misaligned images
imshowpair(fixed(:,:,30), moving(:,:,30),'Scaling','joint');

%% Get a configuration suitable for registering images from different
% sensors.
[optimizer, metric] = imregconfig('monomodal');

%% Find geometric transformation that maps moving to fixed.
tform = imregtform(moving, fixed, 'affine', optimizer, metric);

%% Use imwarp to apply tform to moving so that it aligns with fixed.
% Preserve world limits and resolution of the fixed image when forming
% the transformed image using the 'OutputView' Name/Value pair.
movingRegistered = imwarp(moving,tform,'OutputView',imref3d(size(fixed)));

%% View registered images
figure
imshowpair(fixed(:,:,30), movingRegistered(:,:,30),'Scaling','joint');

tform = {};

movingRegistered = zeros(size(magnitudeData));
angleRegistered = zeros(size(angleData));

%%
for i = 1:size(magnitudeData,4)
    
    moving = magnitudeData(:,:,:,i);
    
    % Find geometric transformation that maps moving to fixed.
    tform_cell{i} = imregtform(moving, fixed, 'rigid', optimizer, metric);

    % Use imwarp to apply tform to moving so that it aligns with fixed.
    % Preserve world limits and resolution of the fixed image when forming
    % the transformed image using the 'OutputView' Name/Value pair.
    movingRegistered(:,:,:,i) = imwarp(moving,tform_cell{i},'OutputView',imref3d(size(fixed)));
    angleRegistered(:,:,:,i) = imwarp(angleData(:,:,:,i),tform_cell{i},'OutputView',imref3d(size(fixed)));
    
end

%%

tform_mat = zeros(4,4,260);

for i = 1:260
tform_mat(:,:,i) = tform_cell{i}.T;
end

figure()
for i = 1:4
    for j = 1:4
   
        plot(squeeze(tform_mat(i,j,:)));
        hold on;
    
    end
end