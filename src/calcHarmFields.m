function [harmfields,mask1,phas,uphas,fl] = calcHarmFields(angleData,mask0,TE,B0,GYRO,vsz)

% Convert from Int16 to phase
phas = double(angleData) ./ 2048 * pi - pi;

% Do the laplacian unwrapping and rescale to units of field
uphas = unwrapLaplacian(phas,mask0,vsz);

for i = 1:size(uphas,4)
    uphas(:,:,:,i) = uphas(:,:,:,i) ./ (B0 * GYRO * TE);
end

% Background field removal
[fl, mask1] = resharp(uphas, mask0, vsz,9:-2*max(vsz):2*max(vsz), 0.05);

% get the harmonic phase evolution
harmfields = uphas - fl;

end

