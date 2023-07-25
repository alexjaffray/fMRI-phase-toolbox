function plotSphericalHarmonics(shmat,time)

PS = PLOT_STANDARDS();

C = distinguishable_colors(16,'w');

% Plot Each Coefficient of Varying Order


figure()
tl = tiledlayout(4,1,'TileSpacing','compact');

% Order 0
nexttile

yyaxis left
ylabel({'Field Fluctuation'; '[ppm]'});
yyaxis right
h1 = plot(time,shmat(:,1),'LineStyle', '-', 'LineWidth', 2);
ax = gca;
ax.ColorOrder = C(1,:);
set(ax,'FontWeight','bold','Box','on','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','YGrid','off','LineWidth',1.5,'FontSize',16);
set(ax.XAxis,'TickDirection','in');
set(ax.YAxis,'TickDirection','out');

ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(ax.YAxis(1),'TickValues',[]);
set(ax,'XTickLabel',[]);

d1 = {'1'};

for i = 1:length(h1)
   
    h1(i).DisplayName = d1{i};
    
end

% Order 1

nexttile

yyaxis left
ylabel({'1^{st} Order'; '[a.u]'});
yyaxis right
h2 = plot(time,shmat(:,2:4),'LineStyle', '-', 'LineWidth', 2);
ax = gca;
ax.ColorOrder = C(2:4,:);
set(ax,'FontWeight','bold','Box','on','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','YGrid','off','LineWidth',1.5,'FontSize',16);
set(ax.XAxis,'TickDirection','in');
set(ax.YAxis,'TickDirection','out');

ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(ax.YAxis(1),'TickValues',[]);
set(ax,'XTickLabel',[]);

d2 = {'x','y','z'};

for i = 1:length(h2)
   
    h2(i).DisplayName = d2{i};
    
end

% Order 2

nexttile

yyaxis left
ylabel({'2^{nd} Order'; '[a.u]'});
yyaxis right
h3 = plot(time,shmat(:,5:9),'LineStyle', '-', 'LineWidth', 2);
ax = gca;
ax.ColorOrder = C(5:9,:);
set(ax,'FontWeight','bold','Box','on','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','YGrid','off','LineWidth',1.5,'FontSize',16);
set(ax.XAxis,'TickDirection','in');
set(ax.YAxis,'TickDirection','out');

ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(ax.YAxis(1),'TickValues',[]);
set(ax,'XTickLabel',[]);

d3 = {'xy','yz','xz','x^2-y^2','2z^2-x^2-y^2'};

for i = 1:length(h3)
   
    h3(i).DisplayName = d3{i};
    
end

% Order 3

nexttile


yyaxis left
ylabel({'3^{rd} Order'; '[a.u]'});
xlabel('Scan Time [s]');
yyaxis right
h4 = plot(time,shmat(:,10:16),'LineStyle', '-', 'LineWidth', 2);
ax = gca;
ax.ColorOrder = C(10:16,:);
set(ax,'FontWeight','bold','Box','on','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','YGrid','off','LineWidth',1.5,'FontSize',16);
set(ax.XAxis,'TickDirection','in');
set(ax.YAxis,'TickDirection','out');

ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(ax.YAxis(1),'TickValues',[]);

d4 = {'xyz','zx^2-zy^2','3yx^2-y^3','(5z^2-(x^2+y^2+z^2))y','(5z^2-(x^2+y^2+z^2))x','5z^3-3(x^2+y^2+z^2)z','x^3-3xy^2'};
% lg4.Layout.Tile = 8;

for i = 1:length(h4)
   
    h4(i).DisplayName = d4{i};
    
end

lg = legend([h1;h2;h3;h4],'Orientation','Horizontal','FontSize',18,'NumColumns',4);
lg.Layout.Tile = "South";

end

