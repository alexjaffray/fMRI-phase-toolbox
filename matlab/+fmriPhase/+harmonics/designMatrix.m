function basis = designMatrix(mask, voxelSize, order)
%DESIGNMATRIX Build real solid-harmonic-like polynomial basis terms.

arguments
    mask
    voxelSize (1,3) double {mustBePositive} = [1 1 1]
    order (1,1) double {mustBeInteger, mustBePositive} = 3
end

if order > 3
    error('fmriPhase:harmonics:UnsupportedOrder', ...
        'Only harmonic orders up to 3 are implemented.');
end

gridSize = size(mask);
xRange = ((1:gridSize(2)) - (gridSize(2) + 1) / 2) * voxelSize(2);
yRange = ((1:gridSize(1)) - (gridSize(1) + 1) / 2) * voxelSize(1);
zRange = ((1:gridSize(3)) - (gridSize(3) + 1) / 2) * voxelSize(3);
[x, y, z] = meshgrid(xRange, yRange, zRange);

x = x(mask);
y = y(mask);
z = z(mask);
one = ones(size(x));

terms = {one, x, y, z};

if order >= 2
    terms = [terms, {x .* y, y .* z, x .* z, x.^2 - y.^2, ...
        2 * z.^2 - x.^2 - y.^2}];
end

if order >= 3
    radiusSquared = x.^2 + y.^2 + z.^2;
    terms = [terms, {x .* y .* z, z .* x.^2 - z .* y.^2, ...
        3 * y .* x.^2 - y.^3, (5 * z.^2 - radiusSquared) .* y, ...
        (5 * z.^2 - radiusSquared) .* x, ...
        5 * z.^3 - 3 * radiusSquared .* z, x.^3 - 3 * x .* y.^2}];
end

basis = cat(2, terms{:});

end
