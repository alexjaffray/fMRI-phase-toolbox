function [respvolTotal,timeVector, rawTotal] = getRespCompSliceTR(inputImage,s,TR,interpolationFactor)

imSize = s(1:3);

imIndex = 0:19:38;

respvolTotal = zeros(size(inputImage));
rawTotal = zeros(19,s(4));

for ii = 1:19

    inputSG = squeeze(inputImage(:,:,imIndex + ii,:));
    imSizeSG = size(inputSG,1:3);
    
    % Decompose the Image using the SVD
    p = reshape(inputSG,prod(size(inputSG,1:3)),s(4));
    [U,S,V] = svd(double(p),"econ");

    % Take the first 5 components of the SVD
    % first component
    componentVector = 1;
    [comp1, raw1] = recomposeSVD(U,S,V,componentVector,imSizeSG);
    zero_ord_vol = comp1;

    % second component
    componentVector = 2;
    [comp2, raw2] = recomposeSVD(U,S,V,componentVector,imSizeSG);

    % third component
    componentVector = 3;
    [comp3, raw3] = recomposeSVD(U,S,V,componentVector,imSizeSG);

    % fourth component
    componentVector = 4;
    [comp4, raw4] = recomposeSVD(U,S,V,componentVector,imSizeSG);

    % fifth component
    componentVector = 5;
    [comp5, raw5] = recomposeSVD(U,S,V,componentVector,imSizeSG);

    % Interpolate whole volume using zero-filling along the time dimension to help visualize things
    if interpolationFactor <= 1
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
        interpAbs = interpft(inputSG,s(4)*interpolationFactor,4);

        % interpolate the 5 SVD components in time
        interpTime1 = interpft(comp1,s(4)*interpolationFactor,4);
        interpTime2 = interpft(comp2,s(4)*interpolationFactor,4);
        interpTime3 = interpft(comp3,s(4)*interpolationFactor,4);
        interpTime4 = interpft(comp4,s(4)*interpolationFactor,4);
        interpTime5 = interpft(comp5,s(4)*interpolationFactor,4);

    else

        interpAbs = inputSG;
        interpTime1 = comp1;
        interpTime2 = comp2;
        interpTime3 = comp3;
        interpTime4 = comp4;
        interpTime5 = comp5;
    end

    % Create the time vector from the interpolated data
    timeVector = 0:TR/interpolationFactor:s(4)*TR;
    timeVector = timeVector(1:s(4)*interpolationFactor);

    d1 = raw1(1,:);
    d2 = raw2(2,:);
    d3 = raw3(3,:);
    d4 = raw4(4,:);
    d5 = raw5(5,:);

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
    
    respvolTotal(:,:,imIndex + ii,:) = respvol;
    rawTotal(ii ,:) = respcomp;

end

