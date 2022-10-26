% Vb_europa_zoom
% Used to generate the plots shown in Figure 2 in
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
fontsize = 8;
linewidth = 2;

%% Zoom (Vb only)
fns = {'..\..\Freezing Simulations\FREZCHEM\v15.1\Cl_Na_Mg_SO4_seawater\Cl_Na_Mg_SO4_seawater',...
    '..\..\Freezing Simulations\FREZCHEM\v15.1\SO4_Mg_Na_Cl\SO4_Mg_Na_Cl'};
figure
for m = 1:length(fns)
    fn = fns{m};
    [T_liq,Sb_liq] = liquidus_FREZCHEM(fn);
    Teut = min(T_liq);
    Smax = 100;
    S = linspace(0,Smax,1000); % ppt
    T = linspace(100,Tmelt(101325),1e3)-273.15; % C
    [P, Tmat, Smat, Vi_V, Vb_V, Vs_V] = volume_fraction_FREZCHEM(T,S,fn);
    
    Vb_V(Tmat<Teut) = NaN;
    Vs_V(Tmat>Teut & Vs_V==0) = NaN;
    
    subplot(1,2,m)  

if m ==1
    ax1 = gca;
    ax1.FontSize = fontsize;
    ax1.TickLabelInterpreter = 'latex';
    yyaxis(ax1,'left')
    h1 = pcolor(ax1,Smat,Tmat,Vb_V);
    h1.EdgeColor = 'none';
    h1.FaceColor = 'interp';
    axis tight
    hold on
    plot(Sb_liq(Sb_liq<=100),T_liq(Sb_liq<=100),'k','linewidth',linewidth)
    [c,h] = contour(Smat,Tmat,Vb_V,[0.06 0.06],'k','linewidth',1);
    clabel(c,h,'fontsize',fontsize,'interpreter','latex')
    ax1.YDir = 'reverse';
    ylim([Teut 0])
    xlabel('Bulk Salinity, $S$ (ppt)','fontsize',fontsize,'interpreter','latex')
    ylabel('Temperature, $T$ ($^{\circ}$C)','fontsize',fontsize,'interpreter','latex')
    title('Chloride-Dominated','fontsize',fontsize,'interpreter','latex')
    colormap(ax1,brewermap([],'Blues'))
    cb1 = colorbar(ax1,'FontSize',fontsize,...
    'TickLabelInterpreter','latex');
    set(get(cb1,'title'),'string','$V_{b}/V$','interpreter','latex',...
    'fontsize',fontsize);
    yyaxis(ax1,'right')
    ylim([1-(max(T)-Teut/(max(T)-min(T))) 1])
    ax1.YDir = 'reverse';
    ax1.YAxis(1).Color = 'k';
    ax1.YAxis(2).Color = 'k';
    ylabel('Conductive Depth (Normalized)','fontsize',fontsize,'interpreter','latex')
else
    ax2 = gca;
    ax2.FontSize = fontsize;
    ax2.TickLabelInterpreter = 'latex';
    yyaxis(ax2,'left')
    h2 = pcolor(ax2,Smat,Tmat,Vb_V);
    h2.EdgeColor = 'none';
    h2.FaceColor = 'interp';
    axis tight
    hold on
    plot(Sb_liq(Sb_liq<=100),T_liq(Sb_liq<=100),'k','linewidth',linewidth)
    [c,h] = contour(Smat,Tmat,Vb_V,[0.06 0.06],'k','linewidth',1);
    clabel(c,h,'fontsize',fontsize,'interpreter','latex')
    ax2.YDir = 'reverse';
    ylim([Teut 0])
    ylabel('Temperature, $T$ ($^{\circ}$C)','fontsize',fontsize,'interpreter','latex')
    xlabel('Bulk Salinity, $S$ (ppt)','fontsize',fontsize,'interpreter','latex')
    title('Sulfate-Dominated','fontsize',fontsize,'interpreter','latex')
    colormap(ax2,brewermap([],'Blues'))
    cb2 = colorbar(ax2,'FontSize',fontsize,...
    'TickLabelInterpreter','latex');
    set(get(cb2,'title'),'string','$V_{b}/V$','interpreter','latex',...
    'fontsize',fontsize);
    yyaxis(ax2,'right')
    ylim([1-(max(T)-Teut/(max(T)-min(T))) 1])
    ax2.YDir = 'reverse';
    ax2.YAxis(1).Color = 'k';
    ax2.YAxis(2).Color = 'k';
    ylabel('Conductive Depth (Normalized)','fontsize',fontsize,'interpreter','latex')
end

end

set(ax1, 'Layer', 'top')
set(ax2, 'Layer', 'top')

%% Figure Formatting
%  [left bottom width height]

% left margin
margin = 0.10;
ax1.Position(1) = margin/2;

% height
htot = 1-2*margin;
ax1.Position(4) = htot;
ax2.Position(4) = htot;

% width
wtot =  1-margin*3.75;
ax1.Position(3) = wtot/2;
ax2.Position(3) = wtot/2;

% left position
ax2.Position(1) = ax1.Position(1)+ax1.Position(3)+2*margin;

% % bottom
ax1.Position(2) = 1.25*margin;
ax2.Position(2) = 1.25*margin;

f = gcf;
f.Color = 'w';
f.Units = 'centimeters';
width = 190/10;
height = 70/10;
f.Position(3:4) = [width height];

