%% Do projection

% define a regular grid covering the domain of the data (arb scale)
[xxmat,yymat,zzmat] = meshgrid(-47.5:1:47.5,-47.5:1:47.5,-28:1:28);
order0 = ones(96,96);

xx0 = xxmat(:,:,1);
xx = xx0(:);
yy0 = yymat(:,:,1);
yy = yy0(:);

slicetrack = repmat(1:57,1,260)';
se = strel('cube',7);

mask = imerode(mask0,se);

testharmfields = reshape(respvol,96,96,[]);
% testharmfields = reshape(harmfields,96,96,[]);
C = zeros(numel(slicetrack),16);

elvec = 1:3:96;

dmat = [];

tic
for sliceidx = 1:numel(slicetrack)

xx = xx0(mask(elvec,elvec,slicetrack(sliceidx)));
yy = yy0(mask(elvec,elvec,slicetrack(sliceidx)));
zz = zzmat(elvec,elvec,slicetrack(sliceidx));
zz = zz(mask(elvec,elvec,slicetrack(sliceidx)));
onevec = order0(mask(elvec,elvec,slicetrack(sliceidx)));

dmat = [onevec xx yy zz xx.*yy yy.*zz xx.*zz xx.^2 - yy.^2 2*zz.^2 - xx.^2 - yy.^2 xx.*yy.*zz zz.*xx.^2 - zz.*yy.^2 3*yy.*xx.^2 - yy.^3 (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*yy (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*xx 5*zz.^3-3*(xx.^2 + yy.^2 + zz.^2).*zz xx.^3 - 3*xx.*yy.^2];

C(sliceidx,:) = basisExpansion(testharmfields(elvec,elvec,sliceidx),dmat,mask(elvec,elvec,slicetrack(sliceidx)));

end
toc

%% Ignore slice timing -> seems more promising :) 
C2 = [];
testharmfields2 = respvol;
xx = xxmat(mask);
yy = yymat(mask);
zz = zzmat(mask);
order = ones(96,96,57);

onevec = order(mask);

dmat = [onevec xx yy zz xx.*yy yy.*zz xx.*zz xx.^2 - yy.^2 2*zz.^2 - xx.^2 - yy.^2 xx.*yy.*zz zz.*xx.^2 - zz.*yy.^2 3*yy.*xx.^2 - yy.^3 (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*yy (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*xx 5*zz.^3-3*(xx.^2 + yy.^2 + zz.^2).*zz xx.^3 - 3*xx.*yy.^2];

tic
for volidx = 1:260

C2(volidx,:) = basisExpansion(testharmfields2(:,:,:,volidx),dmat,mask);

end

toc

%% Plot Each Coefficient of Varying Order
figure()
subplot(4,1,1);
plot(C2(:,1));
legend('1');
subplot(4,1,2);
plot(C2(:,2:4));
legend('x','y','z');
subplot(4,1,3);
plot(C2(:,5:9));
legend('xy','yz','xz','x^2-y^2','2z^2-x^2-y^2');
subplot(4,1,4);
plot(C2(:,10:16));
legend('xyz','zx^2-zy^2','3yx^2-y^3','(5z^2-(x^2+y^2+z^2))y','(5z^2-(x^2+y^2+z^2))x','5z^3-3(x^2+y^2+z^2)z','x^3-3xy^2');

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

outputField = reshape(C2*dmat2',size(C2,1),20,20,12);

%% Do Plotting

d = {};
locs = {};

for i = 1:16
    
    [d{i}, locs{i}] = findpeaks(squeeze(C(:,i)));
    
end