%% PLOTTING FOR PAPER v2 Proc
close all;
dataList = ["data/christina.mat","data/M13.mat","data/M20.mat","data/M20Post.mat","data/M35.mat","data/M35Post.mat","data/M24.mat","data/M37.mat","data/M40.mat","data/M40Post.mat"];
phaseVec = [-1 -1 1 1 -1 -1 1 1 -1 -1];
scale = zeros(10,260);

PS = PLOT_STANDARDS();

%%
figure();
for i = 1:length(dataList)
    
    tmp = load(dataList(i));
    subplot(length(dataList),1,i);
    
    plot(tmp.volTime,tmp.respcompVOL/max(tmp.respcompVOL),'LineStyle', '-', 'LineWidth', 2, 'Color',PS.Blue4);
    hold on;
    plot(tmp.physLogTime,tmp.physLogResp, 'LineStyle', '-.', 'LineWidth', 2, 'Color', PS.MyRed);
    yyaxis right
    ax = gca;
    set(ax,'FontWeight','bold','Box','on','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','YGrid','off','LineWidth',1.5,'ylim',[-1.2 1.2],'FontSize',18);
    set(ax.XAxis(1),'TickDirection','in');
    set(ax.YAxis(1),'TickDirection','out');
    
    yyaxis left
    ylabel(strjoin({'#',num2str(i)}));
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    set(ax.YAxis(1),'TickValues',[]);

    if i ~= 10
        
        set(ax,'XTickLabel',[]);
        
            
    end

    
    if i == 1
        
        legend('Proposed Method (Volume TR)','Breathing Belt','Location','northeast');
        
    end
    
    if i == 10
        
        xlabel('Scan Time [s]');
    
    end
    
    
        
end
%%

peakCell = {};
figure();

for i = 1:length(dataList)
   
    tmp = load(dataList(i));
    tmpResp = tmp.respcompVOL/max(tmp.respcompVOL);

    subplot(10,1,i);
    [f,t] = findpeaks(tmp.physLogResp,tmp.physLogTime,'MinPeakProminence',0.15)
    findpeaks(tmp.physLogResp,tmp.physLogTime,'MinPeakProminence',0.15)
    hold on
    [f2,t2] = findpeaks(-1*tmpResp,tmp.volTime,'MinPeakProminence',0.15)
    findpeaks(tmpResp,tmp.volTime,'MinPeakProminence',0.15)
    
    
    peakCell{i,1} = f;
    peakCell{i,2} = t';
    peakCell2{i,1} = f2;
    peakCell2{i,2} = t2';
end

%%

figure();
for i = 1:length(dataList)
    
    tmp = load(dataList(i));
    subplot(length(dataList),1,i);
    
    plot(tmp.sliceTime,phaseVec(i)*tmp.respcompSL/max(tmp.respcompSL),'LineStyle','-','LineWidth',2,'Color',PS.Green4);
    hold on;
    plot(tmp.physLogTime,tmp.physLogResp, 'LineStyle', '-.', 'LineWidth', 2, 'Color', PS.MyRed);
    yyaxis right
    ax = gca;
    set(ax,'FontWeight','bold','Box','on','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','YGrid','off','LineWidth',1.5,'ylim',[-1.2 1.2],'FontSize',18);
    set(ax.XAxis(1),'TickDirection','in');
    set(ax.YAxis(1),'TickDirection','out');

    yyaxis left
    ylabel(strjoin({'#',num2str(i)}));
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    set(ax.YAxis(1),'TickValues',[]);

    if i ~= 10
        
        set(ax,'XTickLabel',[]);
       
    end
    
    if i == 1
        
        legend('Proposed Method (Slice TR)','Breathing Belt','Location','northeast');

    end
    
    if i == 10
        
        xlabel('Scan Time [s]');
    
    end
end

%% Set the TR
TR = 1.15
interpolationFactor = 1
n_timepoints = 260
n_exc = 19


%% Plot Spherical Harmonics for all datasets
close all;
for i = 1
    
    tmp = load(dataList(i));
    % Plot Each Coefficient of Varying Order (volume TR)
    plotSphericalHarmonics(tmp.volTR_coeffs,TR/interpolationFactor*(1:n_timepoints));
    % Plot Each Coefficient of Varying Order (Slice TR, Filtered Notch + LP)
    plotSphericalHarmonics(tmp.sliceTR_coeffs,TR/interpolationFactor/n_exc*(1:(n_timepoints*n_exc)));
    
    
end

%% Plot depiction of harmonics
% Generate small grid of x y and z positions to visualize the spherical harmonic evolution
tmp = load(dataList(1));

[xxv,yyv,zzv] = meshgrid(-47.5:.5:47.5,-47.5:.5:47.5,-28:0.475:28);

xx = xxv(:);
yy = yyv(:);
zz = zzv(:);

dmat2 = [ones(size(xx)) xx yy zz xx.*yy yy.*zz xx.*zz xx.^2 - yy.^2 2*zz.^2 - xx.^2 - yy.^2 xx.*yy.*zz zz.*xx.^2 - zz.*yy.^2 3*yy.*xx.^2 - yy.^3 (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*yy (5*zz.^2 - (xx.^2 + yy.^2 + zz.^2)).*xx 5*zz.^3-3*(xx.^2 + yy.^2 + zz.^2).*zz xx.^3 - 3*xx.*yy.^2];

n_interp_pts = 1000;
fieldVec = tmp.volTR_coeffs*dmat2';
outputField = reshape(fieldVec,size(tmp.volTR_coeffs,1),191,191,118);
smallMask = imresize3(mask0,[191,191,118]);

% select timepoints

demofig = figure();
timePoints = 104:111;
selslice = 70;
xco = 40:150;
yco = 28:160;

for i = 1:length(timePoints)
   
    subplot(3,4,i);
    
    intermIm = squeeze(outputField(timePoints(i),yco,xco,selslice)).*squeeze(smallMask(yco,xco,selslice));
    intermMask = 1 - smallMask(yco,xco,selslice);
    intermIm2 = imoverlay(intermIm,intermMask,'black');
    
    imagesc(squeeze(outputField(timePoints(i),yco,xco,selslice)).*squeeze(smallMask(yco,xco,selslice)),[-.015 0.015]);
    colormap('turbo');
    ax = gca;
    set(ax,'FontWeight','bold','Box','on','XMinorTick','off','YMinorTick','off','YGrid','off','LineWidth',1.5, 'FontSize',18);
    set(ax.XAxis,'TickDirection','in');
    set(ax.YAxis,'TickDirection','out');
    set(ax,'XTick',[]);
    set(ax,'YTick',[]);
    set(ax,'YTickLabel',[]);
    set(ax,'XTickLabel',[]);
    title(strjoin({'t = ',num2str(timePoints(i)*1.15),'s'}));
end

subplot(3,4,[9,10,11,12]);
plot(tmp.volTime, tmp.respcompVOL,'k','LineWidth',2,'LineStyle','-.');
hold on
plot(tmp.volTime(timePoints),tmp.respcompVOL(timePoints),'r','LineWidth',2,'LineStyle','-','Marker','o');
ax = gca;
set(ax,'FontWeight','bold','Box','on','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','YGrid','off','LineWidth',1.5, 'FontSize',18);
set(ax.XAxis,'TickDirection','in');
set(ax.YAxis,'TickDirection','out');
xlabel('Scan Time [s]');
ylabel('B_0 Fluctuation [ppm]');
xline(tmp.volTime(timePoints(1)),'LineStyle','-.','LineWidth',1.5);
xline(tmp.volTime(timePoints(end)),'LineStyle','-.','LineWidth',1.5);
legend('0th Order B_0 Fluctuation','Sample Times (Shown)');

demoh = axes(demofig,'visible','off');

c = colorbar(demoh,'Position',[0.105 0.42 0.01 0.5],'FontSize',18,'FontWeight','bold');
c.Label.String = "B_0 [ppm]";
colormap(c,'turbo');
caxis(demoh,[-0.015 0.015]);

% %% Write Video
% 
% inField = interp1(tmp.volTime,plat,linspace(tmp.volTime(1),tmp.volTime(end),n_interp_pts),'cubic');
% outputField = reshape(inField,n_interp_pts,191,191,118);
% 
% V = VideoWriter('myFile2.avi');
% V.VideoCompressionMethod
% open(V)
% 
% imagesc(squeeze(outputField(100,:,:,selslice)).*squeeze(smallMask(:,:,selslice)),[-0.02,0.02])
% 
% %% 
% colormap(turbo)
% axis tight manual
% set(gca,'nextplot','replacechildren');
% for k = 1:(size(outputField,1))
%     imagesc(squeeze(outputField(k,:,:,selslice)).*squeeze(smallMask(:,:,selslice)))
%     frame = getframe(gcf);
%     writeVideo(V,frame);
%     
%     k
%     
% end
% 
% close(V);

%% Get Group Metrics

%% Plot Spherical Harmonics for all datasets

metExc = zeros(1,10);

close all;
for i = 1:length(dataList)
    
    tmp = load(dataList(i));
    % Plot Each Coefficient of Varying Order (volume TR)
    % plotSphericalHarmonics(tmp.volTR_coeffs,TR/interpolationFactor*(1:n_timepoints));
    % Plot Each Coefficient of Varying Order (Slice TR, Filtered Notch + LP)
    % plotSphericalHarmonics(tmp.sliceTR_coeffs,TR/interpolationFactor/n_exc*(1:(n_timepoints*n_exc)));
    
    fieldVec = tmp.volTR_coeffs*dmat2';
    outputField = reshape(fieldVec,size(tmp.volTR_coeffs,1),191,191,118);
    maxVal  = 0;
    minVal = 0;
    size(outputField)
    for ll = 1:260
        
        t23 = max(squeeze(outputField(ll,:,:,:)).*smallMask,[],'all');
        t21 = min(squeeze(outputField(ll,:,:,:)).*smallMask,[],'all');
        
        if t23 > maxVal
            maxVal = t23;
        end
        if t21 < minVal
            minVal = t21;
        end
        
    end
    
    metExc(i) = maxVal - minVal;
    
    
end


