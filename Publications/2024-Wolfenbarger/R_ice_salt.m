% R_ice_salt.m 
% Used to generate the plot shown in Figure 1a in
%
% Wolfenbarger, N. S., Blankenship, D. D., Young, D. A., Scanlan, K. M.,
% Chivers, C. J., Findlay, D., Steinbruegge, G. B., Chan, K., Grima, C.,
% Soderlund, K. M., & Schroeder, D. M. (XXXX). Radar Characterization of
% Salt Layers in Europa’s Ice Shell as a Window into Critical Ice-Ocean
% Exchange Processes.

clear all; close all; clc
addpath('..\..\PHREEQC')
addpath('..\..\Ice Ih Properties')
addpath('..\..\Water Properties')

%% Defaults
fontsize = 10;
linewidth = 2;

%% Plot
figure
fn = '..\..\Freezing Simulations\PHREEQC\frezchem_ColdChem\MgSO4_1ppt.pqo';

[~,T_liq,Sb_liq] = liquidus_PHREEQC(fn);
Seut = max(Sb_liq);
S = linspace(0,Seut,1000);
T = linspace(50,160)-273.15;


[P, Tmat, Smat, Vi_V, Vb_V, Vs_V] = volume_fraction_PHREEQC(T,S,fn);
Vs_V_mean = mean(Vs_V,1)';
Vs_V_min = min(Vs_V,[],1)';
Vs_V_max = max(Vs_V,[],1)';

patch([S'; flip(S')], [Vs_V_min; flip(Vs_V_max)], [1 1 1]*0.8, 'EdgeColor','none')
hold on
plot(S,Vs_V_mean,'k','LineWidth',linewidth);

ax = gca;
axis tight
ax.FontSize = fontsize;
ax.XLabel.String = 'Bulk Ice Salinity, $S$ (ppt)';
ax.XLabel.FontSize = fontsize;
ax.YLabel.String = 'Salt Volume Fraction, $V_{ss}/V$';
ax.YLabel.FontSize = fontsize;
ax.XMinorTick = 'on';
box on


%% Figure Formatting
% Figure Formatting
%  [left bottom width height]
M = 1;
N = 2;

% left
margin = 0.175;
ax.Position(1) = margin/N*N;

% height
htot = 1-(margin/M+margin/2);
ax.Position(4) = htot/M;

% width
wtot = 1-(margin/N+3/4*margin);
ax.Position(3) = wtot;

% bottom
ax.Position(2) = margin/M;

f = gcf;
f.Color = 'w';
f.Units = 'inches';
width = 3;
height = 3;
f.Position(3:4) = [width height];




