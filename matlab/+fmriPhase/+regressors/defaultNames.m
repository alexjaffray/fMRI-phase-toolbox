function names = defaultNames(nHarmonic, nSvd)
%DEFAULTNAMES Return stable, package-safe regressor column names.

arguments
    nHarmonic (1,1) double {mustBeInteger, mustBeNonnegative}
    nSvd (1,1) double {mustBeInteger, mustBeNonnegative}
end

names = strings(1, nHarmonic + nSvd);

for idx = 1:nHarmonic
    names(idx) = sprintf('phase_harmonic_%02d', idx);
end

for idx = 1:nSvd
    names(nHarmonic + idx) = sprintf('phase_svd_%02d', idx);
end

end
