function [f, t] = getROImean(vol,roix,roiy,roiz,td)
%getROImean gets the time evolution of the region of interest inside of vol
%  vol: 4D data [3d spatial] x [1d time] -> [Nx Ny Nz Nt]
%  roi: vector with range of indices [3d spatial] -> [3]
%  f: mean of vol(roi) at each time point
%  t: time points (0:td:td*Nt)

f = squeeze(mean(vol(roix,roiy,roiz,:),[1 2 3]));
t = 0:td:td*size(vol,4);

end

