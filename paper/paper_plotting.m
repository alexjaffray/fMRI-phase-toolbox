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

