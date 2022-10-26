% habitability_metrics_europa.m
% Used to generate the plots shown in Figure 1 in
%
% Wolfenbarger, N. S., Fox-Powell, M. G., Buffo, J. J., Soderlund, K. M., &
% Blankenship, D. D. (submitted). The Habitability of Brine Pockets in
% Europa's Ice Shell. Geophysical Research Letters,
% https://doi.org/10.1002/essoar.10512037.1.


clear all; close all; clc
addpath('..\..\FREZCHEM')
addpath('..\..\Ice Ih Properties')
addpath('..\..\Water Properties')

%% Defaults
colors = brewermap(10,'Set2');
fdir = '..\..\Freezing Simulations\FREZCHEM\v15.1\';
fontsize = 8;

%% Comparison of H2O activity, brine salinity, ionic strength f(Vb/V)
solutions = {'Cl_Na_Mg_SO4_seawater\Cl_Na_Mg_SO4_seawater',...
    'SO4_Mg_Na_Cl\SO4_Mg_Na_Cl',...
    'NaCl\NaCl_1ppt',...
    'MgSO4\MgSO4_1ppt'};

figure
subplot(2,4,1)
plot([0.6 1],[-10.9 -10.9],'k','linewidth',1)
hold on
plot([0.9 0.9],[0 -35],'k','linewidth',1)
subplot(2,4,2)
plot([0 10],[-10.9 -10.9],'k','linewidth',1)
hold on
subplot(2,4,3)
plot([0 300],[-10.9 -10.9],'k','linewidth',1)
hold on

subplot(2,4,5)
semilogy([0.9 0.9],[1e-4 1],'k','linewidth',1)
hold on

for n = 1:length(solutions)
    fn = solutions{n};
    FrOut = read_FrOut([fdir,fn]);
    Tm = Tmelt(FrOut.pressure(1)*1e5)-273.15; % C
    
    ice_H2O = FrOut.ice_H2O;
    ind = find(ice_H2O>0);
    
    T = FrOut.T(ind)-273.15; % C
    br_H2O = FrOut.br_H2O(ind);
    ms_b = FrOut.br_salt(ind);
    
    Sb = 1000*ms_b./(ms_b+br_H2O);
    H2O_activity = FrOut.H2O_activity(ind);
    ion_str = FrOut.ion_str(ind);
    
    % Extend to T = Tm
    T = [Tm; T];
    Sb = [0; Sb];
    H2O_activity = [1; H2O_activity];
    ion_str = [0; ion_str];
    
    if contains(solutions{n},'seawater')
        aH2O = H2O_activity;
        T_aH2O = T;
    end
    
    
    S = 1;
    [~,~,~,~, Vb, ~] = volume_fraction_FREZCHEM(T,S,[fdir,fn]);
    
    subplot(2,4,1)
    p(n) = plot(H2O_activity,T,'linewidth',2,'color',colors(n,:));
    
    subplot(2,4,2)
    plot(ion_str,T,'linewidth',2,'color',colors(n,:))
    
    subplot(2,4,3)
    plot(Sb,T,'linewidth',2,'color',colors(n,:))
    
    subplot(2,4,5)
    semilogy(H2O_activity,Vb,'linewidth',2,'color',colors(n,:));
    
    subplot(2,4,6)
    semilogy(ion_str,Vb,'linewidth',2,'color',colors(n,:))
    hold on

    subplot(2,4,7)
    semilogy(Sb,Vb,'linewidth',2,'color',colors(n,:))
    hold on

    subplot(2,4,8)
    semilogy(T,Vb,'linewidth',2,'color',colors(n,:))
    hold on

    
    Teut(n) = min(T);
    
    if contains(fn,'SO4_Mg_Na_Cl')
        leg{n} = sprintf('Sulfate-\nDominated');
    elseif contains(fn,'Cl_Na_Mg_SO4_seawater')
        leg{n} = sprintf('Chloride-\nDominated');
    elseif contains(fn,'MgSO4_1ppt')
        leg{n} = 'MgSO4';
    elseif contains(fn,'NaCl_1ppt')
        leg{n} = 'NaCl';
    end
end

subplot(2,4,1)
ax1 = gca;
ax1.Color = 'w';
ax1.XTickLabel = [];
ax1.YLabel.String = 'Temperature, $T$ ($^{\circ}$C)';
ax1.YLabel.Interpreter = 'latex';
ax1.XLim = [0.6 1];
ax1.YLim = [0.8*(273.15-100)+100  273.15]-273.15;
ax1.YDir = 'reverse';
ax1.FontSize = fontsize;
ax1.TickLabelInterpreter = 'latex';


subplot(2,4,2)
ax2 = gca;
ax2.Color = 'w';
ax2.XTickLabel = [];
ax2.XLim = [0 10];
ax2.FontSize = fontsize;
ax2.TickLabelInterpreter = 'latex';
ax2.YLim = [0.8*(273.15-100)+100  273.15]-273.15;
ax2.YTickLabel = [];
ax2.YDir = 'reverse';

subplot(2,4,3)
ax3 = gca;
ax3.XTickLabel = [];
ax3.XLim = [0 300];
ax3.FontSize = fontsize;
ax3.TickLabelInterpreter = 'latex';
yyaxis(ax3,'left')
ax3.YLim = [0.8*(273.15-100)+100  273.15]-273.15;
ax3.YTickLabel = [];
ax3.YDir = 'reverse';
yyaxis(ax3,'right')
ax3.YLim = [0.8 1];
ax3.YDir = 'reverse';
ax3.YLabel.String = 'Conductive Depth (Normalized)';
ax3.YLabel.Interpreter = 'latex';
ax3.YAxis(2).Color = 'k';

subplot(2,4,4)
axis off

subplot(2,4,5)
ax4 = gca;
ax4.XLabel.String = 'Water Activity';
ax4.YLabel.String = 'Brine Volume Fraction, $V_b/V$';
ax4.XLabel.Interpreter = 'latex';
ax4.YLabel.Interpreter = 'latex';
ax4.XLim = [0.6 1];
ax4.YLim = [1e-4 1];
ax4.YTick = [1e-4 1e-3 1e-2 1e-1 1];
ax4.YDir = 'reverse';
ax4.FontSize = fontsize;
ax4.TickLabelInterpreter = 'latex';

subplot(2,4,6)
ax5 = gca;
ax5.XLabel.String = {'Ionic Strength','(mol/kg H$_{2}$O)'};
ax5.XLabel.Interpreter = 'latex';
ax5.XLim = [0 10];
ax5.YTickLabel = [];
ax5.YLim = [1e-4 1];
ax5.YDir = 'reverse';
ax5.FontSize = fontsize;
ax5.TickLabelInterpreter = 'latex';

subplot(2,4,7)
ax6 = gca;
ax6.XLabel.String = {'Brine Salinity, $S_b$','(ppt)'};
ax6.XLabel.Interpreter = 'latex';
ax6.XLim = [0 300];
ax6.YTickLabel = [];
ax6.YLim = [1e-4 1];
ax6.YDir = 'reverse';
ax6.FontSize = fontsize;
ax6.TickLabelInterpreter = 'latex';

subplot(2,4,8)
ax7 = gca;
ax7.XLabel.String = {'Temperature, $T$','($^{\circ}$C)'};
ax7.XLabel.Interpreter = 'latex';
ax7.XLim = [-35 0];
ax7.YTickLabel = [];
ax7.YLim = [1e-4 1];
ax7.YDir = 'reverse';
ax7.FontSize = fontsize;
ax7.TickLabelInterpreter = 'latex';

%% Figure Formatting
%  [left bottom width height]
M = 2;
N = 4;

% left
margin = 0.175;
ax1.Position(1) = margin/3;
ax4.Position(1) = margin/3;

% height
htot = 1-(margin/M+margin/2);
ax1.Position(4) = htot/M;
ax2.Position(4) = htot/M;
ax3.Position(4) = htot/M;
ax4.Position(4) = htot/M;
ax5.Position(4) = htot/M;
ax6.Position(4) = htot/M;
ax7.Position(4) = htot/M;


% width
delta_w = margin/16;
wtot = 1-(margin/N+margin/2);
ax1.Position(3) = wtot/N-delta_w;
ax2.Position(3) = wtot/N-delta_w;
ax3.Position(3) = wtot/N-delta_w;
ax4.Position(3) = wtot/N-delta_w;
ax5.Position(3) = wtot/N-delta_w;
ax6.Position(3) = wtot/N-delta_w;
ax7.Position(3) = wtot/N-delta_w;

% left
ax2.Position(1) = ax1.Position(1)+ax1.Position(3)+margin/N/2; 
ax3.Position(1) = ax2.Position(1)+ax2.Position(3)+margin/N/2; 
ax5.Position(1) = ax4.Position(1)+ax4.Position(3)+margin/N/2; 
ax6.Position(1) = ax5.Position(1)+ax5.Position(3)+margin/N/2; 
ax7.Position(1) = ax6.Position(1)+ax6.Position(3)+margin/N/2; 

% bottom
delta_bot = margin/8;
ax7.Position(2) = margin/M+delta_bot;
ax6.Position(2) = margin/M+delta_bot;
ax5.Position(2) = margin/M+delta_bot;
ax4.Position(2) = margin/M+delta_bot;
ax3.Position(2) = ax7.Position(2)+ax7.Position(4)+margin/M/2;
ax2.Position(2) = ax6.Position(2)+ax6.Position(4)+margin/M/2;
ax1.Position(2) = ax5.Position(2)+ax5.Position(4)+margin/M/2;

% labels
annotation('textbox',[ax1.Position(1)+delta_w ax1.Position(2)+ax1.Position(4)-2*delta_bot 0.1 0.1],...
    'string','a','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')
annotation('textbox',[ax2.Position(1)+delta_w ax2.Position(2)+ax2.Position(4)-2*delta_bot 0.1 0.1],...
    'string','b','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')
annotation('textbox',[ax3.Position(1)+delta_w ax3.Position(2)+ax3.Position(4)-2*delta_bot 0.1 0.1],...
    'string','c','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')
annotation('textbox',[ax4.Position(1)+delta_w ax4.Position(2)+ax4.Position(4)-2*delta_bot 0.1 0.1],...
    'string','d','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')
annotation('textbox',[ax5.Position(1)+delta_w ax5.Position(2)+ax5.Position(4)-2*delta_bot 0.1 0.1],...
    'string','e','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')
annotation('textbox',[ax6.Position(1)+delta_w ax6.Position(2)+ax6.Position(4)-2*delta_bot 0.1 0.1],...
    'string','f','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')
annotation('textbox',[ax7.Position(1)+delta_w ax7.Position(2)+ax7.Position(4)-2*delta_bot 0.1 0.1],...
    'string','g','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')


f = gcf;
f.Color = 'w';
f.Units = 'centimeters';
width = 3/3*190/10;
height = 2/4*230/10;
f.Position(3:4) = [width height];

% legend
leg = legend(p,leg,'fontsize',fontsize,'interpreter','latex');
leg.Position(2) = ax4.Position(2)+ax4.Position(4)+margin/M/2+...
    ax4.Position(4)/2-leg.Position(4)/2;
leg.Position(1) = ax7.Position(1)+ax7.Position(3)-leg.Position(3);

