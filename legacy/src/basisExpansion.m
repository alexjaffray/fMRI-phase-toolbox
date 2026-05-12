function C = basisExpansion(slice,dmat,mask)
%BASISEXPANSION fit solid harmonic basis functions to field
%   

% best-fit quadratic curve
C = dmat \ slice(mask);


end

