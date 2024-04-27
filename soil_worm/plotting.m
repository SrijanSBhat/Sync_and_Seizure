clc; close all; clear

load('workspace.mat')

%%
temp = zeros(277);
temp((A == 1)) = 0.6;

h = imagesc(temp);
colormap(flipud(slanCM('hot')));
% set(h, 'EdgeColor', 'none','LineWidth',0.0001);
set(gca,'TickLabelInterpreter','latex','FontSize',17)
set(gca,'LineWidth',0.8)
xlim([1 277])
ylim([1 277])
% colorbar
% clim([-0.8 1])
%% Spatial-temporal evolution
close all
%

i = 5795;
ii = 1000 + i;

figure('Position', [40 40 1000 400])
subplot(1,3,1)
%
scatter(p(ii,c1),c1,'MarkerFaceColor','#01CAFF','MarkerEdgeColor','none')

hold on

scatter(p(ii,c2),c2,'MarkerFaceColor','#5EC611','MarkerEdgeColor','none')
scatter(p(ii,c3),c3,'MarkerFaceColor','#EA8615','MarkerEdgeColor','none')
scatter(p(ii,c4),c4,'MarkerFaceColor','#C689E7','MarkerEdgeColor','none')
scatter(p(ii,c5),c5,'MarkerFaceColor','#FD5E83','MarkerEdgeColor','none')
scatter(p(ii,c6),c6,'MarkerFaceColor','#2D9173','MarkerEdgeColor','none')

set(gca,'TickLabelInterpreter','latex','FontSize',17)
ax = gca;
set(gca,'LineWidth',0.8)
xlabel('$p_{i}$','Interpreter','latex')
ylabel('$i$','Interpreter','latex')
axis tight
set(ax, "Box","on")

subplot(1,3,[2,3])
xline(i,'-','LineWidth',2,'Color','k')
hold on
imagesc(p(1000:end,:)');
set(gca,'TickLabelInterpreter','latex','FontSize',17)
set(gca,'LineWidth',0.8)
colormap(slanCM('YlOrBr'));
a=colorbar('TickLabelInterpreter','latex');
ylabel(a,'$p_i$','FontSize',16,'Rotation',0,'Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')
xticks([1 2000 4000 6000 8000 10000 12000])
xticklabels({'0','100','200','300','400','500','600'})
set(gca, "Box","on")
axis tight

%% movie

writerObj = VideoWriter('test2.avi'); %// initialize the VideoWriter object

open(writerObj) ;
figure('Position', [40 40 1000 425])
for i = 2000:15:14001


    
    subplot(1,3,1)
    %
    scatter(p(c1),c1,'MarkerFaceColor','#03CAFD','MarkerEdgeColor','none')

    hold on

    scatter(p(i,c2),c2,'MarkerFaceColor','#5EC611','MarkerEdgeColor','none')
    scatter(p(i,c3),c3,'MarkerFaceColor','#FD5E83','MarkerEdgeColor','none')
    scatter(p(i,c4),c4,'MarkerFaceColor','#C68AE7','MarkerEdgeColor','none')
    scatter(p(i,c5),c5,'MarkerFaceColor','#EA8416','MarkerEdgeColor','none')
    scatter(p(i,c6),c6,'MarkerFaceColor','#2D9172','MarkerEdgeColor','none')

    set(gca,'TickLabelInterpreter','latex','FontSize',17)
    ax = gca;
    set(gca,'LineWidth',0.8)
    xlabel('$p_{i}$','Interpreter','latex')
    ylabel('$i$','Interpreter','latex')
    axis tight
    set(ax, "Box","on")

    subplot(1,3,[2,3])
    xline(i,'-','LineWidth',3,'Color','k')
    hold on
    imagesc(p(2000:end,:)');
    set(gca,'TickLabelInterpreter','latex','FontSize',17)
    set(gca,'LineWidth',0.8)
    colormap(slanCM('YlOrBr'));
    a=colorbar('TickLabelInterpreter','latex');
    ylabel(a,'$p_i$','FontSize',16,'Rotation',0,'Interpreter','latex')
    xlabel('$t$ (s)','Interpreter','latex')
    xticks([1 2000 4000 6000 8000 10000 12000])
    xticklabels({'0','100','200','300','400','500','600'})
    set(gca, "Box","on")
    axis tight


    F = getframe ;           %// Capture the frame
    writeVideo(writerObj,F)  %// add the frame to the movie



end
close(writerObj);


%% Gephi

EdgeL = adj2gephilab('worm_network',A,l');


%% Hubness

ER_graph = makerandCIJ_und(277,1918);
k_r = mean(sum(ER_graph,2)); sig_r = std(sum(ER_graph,2));

rho = 2*1918/(N*(N-1));

%h_i = (deg - (N-1)*rho)./sqrt((N-1)*rho*(1-rho));
h_i = (deg - k_r)./sig_r;
%% Participation coeff

P_i = participation_coef(A,l');

%% Hubness - P index plot
close all

figure
hold on
scatter(P_i(c1),h_i(c1),'MarkerFaceColor','#03CAFD','MarkerEdgeColor','#2295B6','SizeData',40,'LineWidth',0.85)
scatter(P_i(c2),h_i(c2),'MarkerFaceColor','#5EC611','MarkerEdgeColor','#45910D','SizeData',40,'LineWidth',0.85)
scatter(P_i(c3),h_i(c3),'MarkerFaceColor','#FD5E83','MarkerEdgeColor','#B95268','SizeData',40,'LineWidth',0.85)
scatter(P_i(c4),h_i(c4),'MarkerFaceColor','#C68AE7','MarkerEdgeColor','#8A5CA3','SizeData',40,'LineWidth',0.85)
scatter(P_i(c5),h_i(c5),'MarkerFaceColor','#EA8416','MarkerEdgeColor','#A05E0E','SizeData',40,'LineWidth',0.85)
scatter(P_i(c6),h_i(c6),'MarkerFaceColor','#2D9172','MarkerEdgeColor','#21644F','SizeData',40,'LineWidth',0.85)
ylim([-5 20])
xlim([0 0.9])
set(gca,'TickLabelInterpreter','latex','FontSize',17)
ax = gca;
set(gca,'LineWidth',0.8)
xlabel('Participation co-eff','Interpreter','latex')
ylabel('Hubness','Interpreter','latex')
%axis tight
set(ax, "Box","on")
legend({'Community 1','Community 2','Community 3','Community 4','Community 5','Community 6'},'EdgeColor','#454955','Interpreter','latex')


%% time series of a neuron
figure('Position', [40 40 1000 400])
plot(p(1000:end,124))

set(gca,'TickLabelInterpreter','latex','FontSize',17)
ax = gca;
set(gca,'LineWidth',0.8)
ylabel('$p_{124}$','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')
axis tight
set(ax, "Box","on")
ylim([-2 2])

xticks([1 2000 4000 6000 8000 10000 12000])
xticklabels({'0','100','200','300','400','500','600'})
