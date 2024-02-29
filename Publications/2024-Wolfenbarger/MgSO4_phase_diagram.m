% MgSO4_phase_diagram.m
% Used to generate the plot shown in Figure 2 in
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
fn = '..\..\Freezing Simulations\PHREEQC\frezchem_ColdChem\MgSO4_1ppt.pqo';

figure
[~,T,Sb] = liquidus_PHREEQC(fn);
plot(Sb,T+273.15,'k','LineWidth',linewidth);

ax = gca;
ax.FontSize = fontsize;
ax.XLabel.String = 'Brine Salinity, $S_b$ (ppt)';
ax.XLabel.FontSize = fontsize;
ax.YLabel.String = 'Temperature, $T$ (K)';
ax.YLabel.FontSize = fontsize;
ax.YLim = [min(T) max(T)]+273.15;


%% Figure Formatting
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
