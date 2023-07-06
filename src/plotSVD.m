function getRespComp(interpAbs,interpAbs)

%% Define Ranges in the slice to look at Fluctuation and calculate roi means
selectedSlice = 31;

width = 8;
height = 10;
xmin = 45;
ymin = 63;

range2 = xmin:(xmin+height);
range1 = ymin:(ymin+width);

d0 = getROImean(interpAbs,range2,range1,selectedSlice,TR/interpolationFactor);
d1 = getROImean(interpTime1,range2,range1,selectedSlice,TR/interpolationFactor);
d2 = getROImean(interpTime2,range2,range1,selectedSlice,TR/interpolationFactor);
d3 = getROImean(interpTime3,range2,range1,selectedSlice,TR/interpolationFactor);
d4 = getROImean(interpTime4,range2,range1,selectedSlice,TR/interpolationFactor);
d5 = getROImean(interpTime5,range2,range1,selectedSlice,TR/interpolationFactor);

%% Plot 3 Images Side by Side
figure();
subplot(2,2,1);
imagesc(phas(:,:,selectedSlice,5));
xlabel('x-position [voxels]');
ylabel('y-position [voxels]');
subplot(2,2,2);
imagesc(uphas(:,:,selectedSlice,5));colormap(gray)
xlabel('x-position [voxels]');
ylabel('y-position [voxels]');
subplot(2,2,3);
imagesc(fl(:,:,selectedSlice,5));colormap(gray)
xlabel('x-position [voxels]');
ylabel('y-position [voxels]');
subplot(2,2,4);
imagesc(harmfields(:,:,selectedSlice,5));colormap(gray)
xlabel('x-position [voxels]');
ylabel('y-position [voxels]');

%% Define axesranges for the plotting
dataMin = min(d0,[],'all');
dataMax = max(d0,[],'all');
dataRange = dataMax - dataMin;
dlim = [0 0];

if dataMin < 0
    
    dlim(1) = dataMin - dataRange;
    dlim(2) = dataRange/2;
    
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
%title('Image with ROI');
xlabel('x-position [voxels]');
ylabel('y-position [voxels]');
colormap(gray);
ax = gca;
hold on;

roi = drawrectangle(ax,'Position',[xmin,ymin,width,height],'Color',[1 0 0]);
subplot(1,2,2);
plot(timeVector,d0,'Color','black','LineStyle','-.');
hold on;
plot(timeVector,d1,'Color','red');
hold on;
plot(timeVector,d2,'Color','blue');
hold on;
plot(timeVector,d3,'Color','black');
hold on;
plot(timeVector,d4,'Color','green');
hold on;
plot(timeVector,d5,'Color','magenta');

%title('Field Fluctuations: Original Image and First 5 SVD Coefficients');
xlabel('Time [s]');
ylabel('Mean Field Fluctuation in ROI [ppm]');

axis([timeVector(1) timeVector(end) dlim(1) dlim(2)]);

legend('Original Image','1st Component','2nd Component','3rd Component','4th Component','5th Component', 'Location','SouthEast');

%% Plot the singular values
figure();
set(gcf,'Position',[100 100 1500 1000]);
plot(1:s(4),S*ones(s(4),1));
title('Singular Values');
xlabel('Singular Value Number');
ylabel('Singular Value');
set(gca,'yscale','log');

%% Determine component of svd corresponding to breathing flux
respcomp = [];

p1 = norm(diff(d1),1);
p2 = norm(diff(d2),1);
p3 = norm(diff(d3),1);
p4 = norm(diff(d4),1);
p5 = norm(diff(d5),1);

respvol = [];

normCheck = 0;
if p1 > normCheck
    respcomp = d1;
    normCheck = p1;
    respvol = interpTime1;
    disp("resp comp = 1")
end
if p2 > normCheck
    respcomp = d2;
    normCheck = p2;
    respvol = interpTime2;
    disp("resp comp = 2")
end
if p3 > normCheck
    respcomp = d3;
    normCheck = p3;
    respvol = interpTime3;
    disp("resp comp = 3")
end
if p4 > normCheck
    respcomp = d4;
    normCheck = p4;
    respvol = interpTime4;
    disp("resp comp = 4")
end
if p5 > normCheck
    respcomp = d5;
    normCheck = p5;
    respvol = interpTime5;
    disp("resp comp = 5")
end

%% Plot the time-course of the 2nd SVD coefficient and identify minima
foundMins = islocalmin(respcomp,'MinProminence',0.0035);

% Plot
figure();
plot(timeVector,respcomp);
hold on;
scatter(timeVector(foundMins),respcomp(foundMins),'o');
title('2nd SVD Coefficient');
xlabel('Time (s)');
ylabel('Fluctuation in X map (ppm)');
legend('Delta X','Minima');
% 
% %% Plot the standard deviation of the qsm values obtained throughout the acquisition
% xRange = std(harmfields,0,4);
% figure(4);
% imagesc(xRange(:,:,25)),colormap('hot');
% colorbar;

end

