% Vb_plots.m 
% Used to generate the plots shown in Figures 7 and 8 in
%
% Wolfenbarger, N. S., Fox-Powell, M. G., Buffo, J. J., Soderlund, K. M., &
% Blankenship, D. D. (2022). Compositional controls on the distribution of
% brine in Europa's ice shell. Journal of Geophysical Research: Planets,
% 127, e2022JE007305. https://doi.org/10.1029/2022JE007305

clear all; close all; clc
addpath('..\..\FREZCHEM')
addpath('..\..\Ice Ih Properties')
addpath('..\..\Water Properties')

%% Defaults
fontsize = 8;
linewidth = 2;

%% Binary
fn = {'..\..\Freezing Simulations\FREZCHEM\v15.1\NaCl\NaCl_1ppt',...
    '..\..\Freezing Simulations\FREZCHEM\v15.1\MgSO4\MgSO4_1ppt'};
sol = {'NaCl','MgSO$_4$'};

for m = 1:length(sol)
    figure(m)
    [T_liq,Sb_liq] = liquidus_FREZCHEM(fn{m});
    Teut = min(T_liq);
    S = linspace(0,max(Sb_liq),100); % ppt
    T = linspace(100,Tmelt(101325),1000)-273.15; % C
    [P, Tmat, Smat, Vi_V, Vb_V, Vs_V] = volume_fraction_FREZCHEM(T,S,fn{m});
    
    Vb_V(Tmat<Teut) = NaN;
    Vs_V(Tmat>Teut & Vs_V==0) = NaN;
    
    ax1 = axes;
    ax1.FontSize = fontsize;
    ax1.TickLabelInterpreter = 'latex';
    yyaxis(ax1,'left')
    h1 = pcolor(ax1,Smat,Tmat,Vb_V);
    h1.EdgeColor = 'none';
    h1.FaceColor = 'interp';
    axis tight
    ax1.YDir = 'reverse';
    xlabel('Bulk Ice Salinity, $S$ (ppt)','fontsize',fontsize,'interpreter','latex')
    ylabel('Temperature, $T$ ($^{\circ}$C)','fontsize',fontsize,'interpreter','latex')
    title(sol{m},'fontsize',fontsize,'interpreter','latex')
    yyaxis(ax1,'right')
    ylim([0 1])
    ax1.YDir = 'reverse';
    ax1.YAxis(1).Color = 'k';
    ax1.YAxis(2).Color = 'k';
    ylabel('Conductive Depth (Normalized)','fontsize',fontsize,'interpreter','latex')
    
    view(2)
    ax2 = axes;
    ax2.FontSize = fontsize;
    ax2.TickLabelInterpreter = 'latex';
    h2 = pcolor(ax2,Smat,Tmat,Vs_V);
    h2.EdgeColor = 'none';
    h2.FaceColor = 'interp';
    axis tight
    hold on
    plot(Sb_liq,T_liq,'k','linewidth',linewidth)
    plot(Sb_liq,Teut*ones(size(T_liq)),'k','linewidth',linewidth/2)
    
    ax2.YDir = 'reverse';
    yyaxis(ax2,'right')
    ylim([0 1])
    ax2.YDir = 'reverse';
    ax2.YAxis(1).Color = 'k';
    ax2.YAxis(2).Color = 'k';
    
    linkaxes([ax1,ax2])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    
    colormap(ax1,brewermap([],'Blues'))
    colormap(ax2,brewermap([],'Reds'))
    set([ax1,ax2],'Position',[.17 .11 .685 .815]);
    cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815],'FontSize',fontsize,...
        'TickLabelInterpreter','latex');
    cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815],'FontSize',fontsize,...
        'TickLabelInterpreter','latex');
    set(get(cb1,'title'),'string','$V_{b}/V$','interpreter','latex',...
        'fontsize',fontsize);
    set(get(cb2,'title'),'string','$V_{ss}/V$','interpreter','latex',...
        'fontsize',fontsize);
    
    set(ax1, 'Layer', 'top')
    set(ax2, 'Layer', 'top')
    
    % Figure Formatting
    % left
    margin = 0.3;
    ax1.Position(1) = margin;
    ax2.Position(1) = margin;
    
    % width
    ax1.Position(3) = 0.4;
    ax2.Position(3) = 0.4;
    cb1.Position(1) = 0.075;
    cb2.Position(1) = 0.85;
    
    f = gcf;
    f.Color = 'w';
    f.Units = 'centimeters';
    width = 95/10;
    height = 115/10;
    f.Position(3:4) = [width height];
end

%% Analog
fn = {'..\..\Freezing Simulations\FREZCHEM\v15.1\Cl_Na_Mg_SO4_seawater\Cl_Na_Mg_SO4_seawater',...
    '..\..\Freezing Simulations\FREZCHEM\v15.1\SO4_Mg_Na_Cl\SO4_Mg_Na_Cl'};
sol = {'Cl-Dominated','SO$_4$-Dominated'};

for m = 1:length(sol)
    figure(m+2)
    [T_liq,Sb_liq] = liquidus_FREZCHEM(fn{m});
    Teut = min(T_liq);
    S = linspace(0,max(Sb_liq),100); % ppt
    T = linspace(100,Tmelt(101325),1e3)-273.15; % C
    [P, Tmat, Smat, Vi_V, Vb_V, Vs_V] = volume_fraction_FREZCHEM(T,S,fn{m});
    
    Vb_V(Tmat<Teut) = NaN;
    Vs_V(Tmat>Teut & Vs_V==0) = NaN;

    ax1 = axes;
    ax1.FontSize = fontsize;
    ax1.TickLabelInterpreter = 'latex';
    yyaxis(ax1,'left')
    h1 = pcolor(ax1,Smat,Tmat,Vb_V);
    h1.EdgeColor = 'none';
    h1.FaceColor = 'interp';
    axis tight
    ax1.YDir = 'reverse';
    xlabel('Bulk Ice Salinity, $S$ (ppt)','fontsize',fontsize,'interpreter','latex')
    ylabel('Temperature, $T$ ($^{\circ}$C)','fontsize',fontsize,'interpreter','latex')
    title(sol{m},'fontsize',fontsize,'interpreter','latex')
    yyaxis(ax1,'right')
    ylim([0 1])
    ax1.YDir = 'reverse';
    ax1.YAxis(1).Color = 'k';
    ax1.YAxis(2).Color = 'k';
    ylabel('Conductive Depth (Normalized)','fontsize',fontsize,'interpreter','latex')
    
    view(2)
    ax2 = axes;
    ax2.FontSize = fontsize;
    ax2.TickLabelInterpreter = 'latex';
    h2 = pcolor(ax2,Smat,Tmat,Vs_V);
    h2.EdgeColor = 'none';
    h2.FaceColor = 'interp';
    h2.FaceAlpha  = '0.75';
    axis tight
    hold on
    plot(Sb_liq,T_liq,'k','linewidth',linewidth)
    plot(Sb_liq,Teut*ones(size(T_liq)),'k','linewidth',linewidth/2)
    
    ax2.YDir = 'reverse';
    yyaxis(ax2,'right')
    ylim([0 1])
    ax2.YDir = 'reverse';
    ax2.YAxis(1).Color = 'k';
    ax2.YAxis(2).Color = 'k';
    
    linkaxes([ax1,ax2])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    
    colormap(ax1,brewermap([],'Blues'))
    colormap(ax2,brewermap([],'Reds'))
    set([ax1,ax2],'Position',[.17 .11 .685 .815]);
    cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815],'FontSize',fontsize,...
        'TickLabelInterpreter','latex');
    cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815],'FontSize',fontsize,...
        'TickLabelInterpreter','latex');
    set(get(cb1,'title'),'string','$V_{b}/V$','interpreter','latex',...
        'fontsize',fontsize);
    set(get(cb2,'title'),'string','$V_{ss}/V$','interpreter','latex',...
        'fontsize',fontsize);
    
    set(ax1, 'Layer', 'top')
    set(ax2, 'Layer', 'top')
    
    % Figure Formatting
    % left
    margin = 0.3;
    ax1.Position(1) = margin;
    ax2.Position(1) = margin;
    
    % width
    ax1.Position(3) = 0.4;
    ax2.Position(3) = 0.4;
    cb1.Position(1) = 0.075;
    cb2.Position(1) = 0.85;
    
    f = gcf;
    f.Color = 'w';
    f.Units = 'centimeters';
    width = 95/10;
    height = 115/10;
    f.Position(3:4) = [width height];
end