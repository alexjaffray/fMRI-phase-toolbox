%% Recompose the image from the Nth component of the SVD
function [recomp, raw] = recomposeSVD(U,S,V,coeffVec,msize)
    mask1 = zeros(size(S));
    mask1(coeffVec,coeffVec) = 1;
    S = S .* mask1;
    raw = S*V';
    recomp = reshape(U*raw,[msize,size(U,2)]);

end