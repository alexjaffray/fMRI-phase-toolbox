%% Recompose the image from the Nth component of the SVD
function recomposed = recomposeSVD(U,S,V,coeffVec,msize)
    mask1 = zeros(size(S));
    mask1(coeffVec,coeffVec) = 1;
    S = S .* mask1;
    recomposed = reshape(U*S*V',[msize,size(U,2)]);

end