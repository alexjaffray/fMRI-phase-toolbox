%% function to project the fMRI phase data (the harmonic data) onto spherical harmonic basis functions defined below:

function timecourse = projectSphericalHarmonics(phas, mask, nsamples)

phas = phas - mean(phas,'all');

phas = phas * 3.0 * 10^-6;

phas = cumtrapz(linspace(0,260,260),phas,4);

s = size(phas);

gridsize = [s(1) s(2) s(3)];
flatsize = prod(gridsize);

fov = [0.25 0.25 0.2];

xrange = linspace(-fov(1)/2,fov(1)/2,gridsize(1));
yrange = linspace(-fov(2)/2,fov(2)/2,gridsize(2));
zrange = linspace(-fov(3)/2,fov(3)/2,gridsize(3));

[x, y, z] = meshgrid(xrange,yrange,zrange);
linmask = reshape(mask,flatsize,1);
approxsamples = randsample(flatsize,nsamples);

approxsamples(~linmask(approxsamples)) = [];

h0 = ones(size(x));

h1 = x;
h2 = y;
h3 = z;

h4 = x .* y;
h5 = z .* y;
h6 = 3 * z.^2 - x.^2 - y.^2 - z.^2;
h7 = x .* z;
h8 = x.^2 - y.^2;

h9 = 3 * y .* x.^2 - y.^3;
h10 = x .* y .* z;
h11 = (5 * z.^2 - x.^2 - y.^2 - z.^2) .* y;
h12 = 5 * z.^3 - 3 * z .* (x.^2 + y.^2 + z.^2);
h13 = (5 * z.^2 - x.^2 - y.^2 - z.^2) .* x;
h14 = x.^2 .* z - y.^2 .* z;
h15 = x.^3 - 3 * x .* y.^2;

basismatrix = [reshape(h0,flatsize,1), reshape(h1,flatsize,1), reshape(h2,flatsize,1), reshape(h3,flatsize,1),reshape(h4,flatsize,1),reshape(h5,flatsize,1),reshape(h6,flatsize,1),reshape(h7,flatsize,1),reshape(h8,flatsize,1),reshape(h9,flatsize,1),reshape(h10,flatsize,1),reshape(h11,flatsize,1),reshape(h12,flatsize,1),reshape(h13,flatsize,1),reshape(h14,flatsize,1),reshape(h15,flatsize,1)];
reducedbasismatrix = basismatrix(approxsamples,:);

Q = (reducedbasismatrix);
P = Q'*Q;

timecourse = zeros(s(4),16);

for i = 1:s(4)
    
    p = reshape(phas(:,:,:,i),flatsize,1);
    preduced = p(approxsamples);
    
    timecourse(i,:) = P\(Q'*preduced);
    
end

end

