%% Do projection

% define a regular grid covering the domain of the data (arb scale)
[xxmat,yymat,zzmat] = meshgrid(-47.5:1:47.5,-47.5:1:47.5,-28:1:28);
order = ones(96,96,57);

n_timepoints = length(respcomp);
mb = 3;
n_exc = size(order,3)/mb;

xx0 = xxmat(:,:,1);
xx = xx0(:);
yy0 = yymat(:,:,1);
yy = yy0(:);

slicetrack = repmat(1:19,1,n_timepoints)';
se = strel('cube',7);

mask = imerode(mask0,se);

testharmfields = respvol;
% testharmfields = reshape(harmfields,96,96,[]);
C = zeros(n_exc,n_timepoints,16);

elvec = 1:96;

dmat = [];

tic
for volidx = 1:n_timepoints
    for ind = 1:n_exc

        zvec = [ind ind+19 ind+38];

        xx = xxmat(mask(elvec,elvec,zvec));
        yy = yymat(mask(elvec,elvec,zvec));
        zz = zzmat(mask(elvec,elvec,zvec));

        onevec = order(mask(elvec,elvec,zvec));

        dmat = [onevec xx yy zz xx.*yy yy.*zz xx.*zz xx.^2 - yy.^2 2*zz.^2 - xx.^2 - yy.^2 xx.*yy.*zz zz.*xx.^2 - zz.*yy.^2 3*yy.*xx.^2 - yy.^3 (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*yy (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*xx 5*zz.^3-3*(xx.^2 + yy.^2 + zz.^2).*zz xx.^3 - 3*xx.*yy.^2];

        C(ind,volidx,:) = basisExpansion(testharmfields(elvec,elvec,zvec,volidx),dmat,mask(elvec,elvec,zvec));
    end
end
toc

C3 = reshape(C,[],16);

%% Filtering
fs = n_exc/TR;
f0 = 1/TR;
fn = fs/2;
freqRatio = f0/fn;

notchWidth = 0.1;

% Compute zeros
notchZeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

% Compute poles
notchPoles = (1-notchWidth) * notchZeros;

b_notch = poly( notchZeros ); %  Get moving average filter coefficients
a_notch = poly( notchPoles ); %  Get autoregressive filter coefficients

% filter signal x
C4 = filter(b_notch,a_notch,C3);
[b_lp,a_lp] = butter(6,freqRatio/2);
C5 = filter(b_lp,a_lp,C4);

%% Ignore slice timing -> seems more promising :) 
C2 = [];
testharmfields2 = respvol;
xx = xxmat(mask);
yy = yymat(mask);
zz = zzmat(mask);

onevec = order(mask);

dmat = [onevec xx yy zz xx.*yy yy.*zz xx.*zz xx.^2 - yy.^2 2*zz.^2 - xx.^2 - yy.^2 xx.*yy.*zz zz.*xx.^2 - zz.*yy.^2 3*yy.*xx.^2 - yy.^3 (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*yy (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*xx 5*zz.^3-3*(xx.^2 + yy.^2 + zz.^2).*zz xx.^3 - 3*xx.*yy.^2];

tic
for volidx = 1:n_timepoints

    C2(volidx,:) = basisExpansion(testharmfields2(:,:,:,volidx),dmat,mask);

end

toc

%% Plot Each Coefficient of Varying Order (volume TR)
plotSphericalHarmonics(C2,TR/interpolationFactor*(1:n_timepoints));

%% Plot Each Coefficient of Varying Order (slice TR)
plotSphericalHarmonics(C3,TR/interpolationFactor/n_exc*(1:(n_timepoints*n_exc)));

%% Plot Each Coefficient of Varying Order (slice TR, Filtered Notch)
plotSphericalHarmonics(C4,TR/interpolationFactor/n_exc*(1:(n_timepoints*n_exc)));

%% Plot Each Coefficient of Varying Order (Slice TR, Filtered Notch + LP)
plotSphericalHarmonics(C5,TR/interpolationFactor/n_exc*(1:(n_timepoints*n_exc)));

%% Generate small grid of x y and z positions to visualize the spherical harmonic evolution

[xxv,yyv,zzv] = meshgrid(-47.5:5:47.5,-47.5:5:47.5,-28:4.75:28);

xx = xxv(:);
yy = yyv(:);
zz = zzv(:);

dmat2 = [ones(size(xx)) xx yy zz xx.*yy yy.*zz xx.*zz xx.^2 - yy.^2 2*zz.^2 - xx.^2 - yy.^2 xx.*yy.*zz zz.*xx.^2 - zz.*yy.^2 3*yy.*xx.^2 - yy.^3 (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*yy (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*xx 5*zz.^3-3*(xx.^2 + yy.^2 + zz.^2).*zz xx.^3 - 3*xx.*yy.^2];

udiag = diag(ones(size(xx,1),1),1);

L = - diag(ones(size(xx,1),1)) + udiag(1:end-1,1:end-1);

sysmat = dmat2*dmat2'+ 0.001*(L*L');

f = pinv(sysmat);

outputField = reshape(C5*dmat2',size(C5,1),20,20,12);

%% Do Plotting

d = {};
locs = {};

for i = 1:16
    
    [d{i}, locs{i}] = findpeaks(squeeze(C(:,i)));
    
end