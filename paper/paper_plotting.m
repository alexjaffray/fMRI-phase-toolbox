%% PLOTTING FOR PAPER

dataList = ["data/christina.mat","data/M13.mat","data/M20.mat","data/M20Post.mat","data/M35.mat","data/M35Post.mat","data/M24.mat","data/M37.mat","data/M40.mat","data/M40Post.mat"];
phaseVec = [-1 1 1 1 -1 -1 1 1 -1 -1];


PS = PLOT_STANDARDS();

figure();
for i = 1:length(dataList)
    
    tmp = load(dataList(i));
    subplot(length(dataList),1,i);
    plot(tmp.fmriTime,phaseVec(i)*tmp.respcomp/max(tmp.respcomp),'LineStyle', '-', 'LineWidth', 2, 'Color',PS.Blue3);
    hold on;
    plot(tmp.physLogTime,tmp.filteredPhysLogTrace, 'LineStyle', '-.', 'LineWidth', 2, 'Color', PS.MyRed);
    set(gca,'Xticklabel',[]) 
    
end

legend('B_0 Fluctuation (Slice TR)','Breathing Belt');

%% PLOTTING FOR PAPER v2 Proc

dataList = ["data/christina.mat","data/M13.mat","data/M20.mat","data/M20Post.mat","data/M35.mat","data/M35Post.mat","data/M24.mat","data/M37.mat","data/M40.mat","data/M40Post.mat"];
phaseVec = [-1 1 1 1 -1 -1 1 1 -1 -1];
scale = zeros(10,260);


PS = PLOT_STANDARDS();
ft = 'Times';
fsz = 12;

figure();
for i = 1:length(dataList)
    
    tmp = load(dataList(i));
    subplot(length(dataList),1,i);
    plot(tmp.volTime,tmp.respcompVOL/max(tmp.respcompVOL),'LineStyle', '-', 'LineWidth', 2, 'Color',PS.Blue4);
    hold on;
    plot(tmp.physLogTime,tmp.physLogResp, 'LineStyle', '-.', 'LineWidth', 2, 'Color', PS.MyRed);
    set(gca,'Xticklabel',[]) 
    set(gca,'FontSize',fsz,'FontName',ft);
    resampPhysLog = interp1(tmp.physLogTime,tmp.physLogResp,tmp.volTime);
    
    
    
    
end

legend('B0 Fluctuation (Volume TR)','Breathing Belt');

figure();
for i = 1:length(dataList)
    
    tmp = load(dataList(i));
    subplot(length(dataList),1,i);
    
    plot(tmp.sliceTime,phaseVec(i)*tmp.respcompSL/max(tmp.respcompSL),'LineStyle','-','LineWidth',2,'Color',PS.Green4);
    hold on;
    plot(tmp.physLogTime,tmp.physLogResp, 'LineStyle', '-.', 'LineWidth', 2, 'Color', PS.MyRed);
    set(gca,'Xticklabel',[]) 
    set(gca,'FontSize',fsz,'FontName',ft);
    
end
legend('B_0 Fluctuation (Slice TR)','Breathing Belt');

%% Set the TR
TR = 1.15
interpolationFactor = 1
n_timepoints = 260
n_exc = 19
%% Plot Each Coefficient of Varying Order (volume TR)
plotSphericalHarmonics(volTR_coeffs,TR/interpolationFactor*(1:n_timepoints));

%% Plot Each Coefficient of Varying Order (Slice TR, Filtered Notch + LP)
plotSphericalHarmonics(sliceTR_coeffs,TR/interpolationFactor/n_exc*(1:(n_timepoints*n_exc)));