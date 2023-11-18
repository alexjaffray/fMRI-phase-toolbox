function [respvol,timeVector, zero_ord_vol] = getRespComp(inputImage,s,TR,interpolationFactor,doPlot)

imSize = s(1:3);

% Decompose the Image using the SVD
p = reshape(inputImage,prod(imSize),s(4));
[U,S,V] = svd(double(p),"econ");

% Take the first 5 components of the SVD
% first component
componentVector = 1;
[comp1, raw1] = recomposeSVD(U,S,V,componentVector,imSize);
zero_ord_vol = comp1;

% second component
componentVector = 2;
[comp2, raw2] = recomposeSVD(U,S,V,componentVector,imSize);

% third component
componentVector = 3;
[comp3, raw3] = recomposeSVD(U,S,V,componentVector,imSize);

% fourth component
componentVector = 4;
[comp4, raw4] = recomposeSVD(U,S,V,componentVector,imSize);

% fifth component
componentVector = 5;
[comp5, raw5] = recomposeSVD(U,S,V,componentVector,imSize);

% Interpolate whole volume using zero-filling along the time dimension to help visualize things
if interpolationFactor < 1
    doInterpolation = false;
else
    doInterpolation = true;
end

interpAbs = [];
interpTime1 = [];
interpTime2 = [];
interpTime3 = [];
interpTime4 = [];
interpTime5 = [];

if doInterpolation
    
    % interpolate the original image in timetimeVector
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

% Create the time vector from the interpolated data
timeVector = 0:TR/interpolationFactor:s(4)*TR;
timeVector = timeVector(1:s(4)*interpolationFactor);

% Define Ranges in the slice to look at Fluctuation and calculate roi means
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

% Determine component of svd corresponding to breathing flux (voxel-free)
respcomp = [];

p1 = norm(diff(raw1(1,:)'),1);
p2 = norm(diff(raw2(2,:)'),1);
p3 = norm(diff(raw3(3,:)'),1);
p4 = norm(diff(raw4(4,:)'),1);
p5 = norm(diff(raw5(5,:)'),1);

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

if doPlot

    % Define axesranges for the plotting
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

    % Plot the first 5 SVD coefficients in the region of interest defined by range2 and range1
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

    % Plot the singular values
    figure();
    set(gcf,'Position',[100 100 1500 1000]);
    plot(1:s(4),S*ones(s(4),1));
    title('Singular Values');
    xlabel('Singular Value Number');
    ylabel('Singular Value');
    set(gca,'yscale','log');

    % Plot the time-course of the 2nd SVD coefficient and identify minima
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

    % Plot the standard deviation of the qsm values obtained throughout the acquisition
    xRange = std(inputImage,0,4);
    figure(4);
    imagesc(xRange(:,:,25)),colormap('hot');
    colorbar;
end

end

