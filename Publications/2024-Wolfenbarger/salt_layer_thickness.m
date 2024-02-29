% salt_layer_thickness.m 
% Used to generate the plot shown in Figure 3 in
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


[P,T,rho_b] = brine_density_PHREEQC(fn);
[P,T,Sb] = liquidus_PHREEQC(fn);

Seut = max(Sb); % ppt
rhoeut = max(rho_b);

mm_MgSO4 = 120.37;
mm_MgSO411H2O = 318.535;
mv_MgSO411H2O = 207.44;

rhos = mm_MgSO411H2O/mv_MgSO411H2O; % g/cm^3

veut = 0.325;
S = linspace(0,Seut);
rho_b = interp1(Sb,rho_b,S);

d = S/1000.*(rho_b/rhos)*mm_MgSO411H2O/mm_MgSO4*(1/veut);

figure
subplot(1,2,1)

d0 = logspace(-3,3,1e3)';
d_salt = d0*S/1000.*(rho_b/rhos)*mm_MgSO411H2O/mm_MgSO4*(1/veut);

x = repmat(S,length(d0),1);
y = repmat(d0,1,length(S));

p = pcolor(x,y,d_salt);
p.LineStyle = 'none';
p.FaceColor = 'interp';
clim([1e-3 1e3])
axis tight

ax = gca;
ax.ColorScale = 'log';
ax.FontSize = fontsize;
ax.XLabel.String = 'Initial Brine Salinity, $S_{b,0}$ (ppt)';
ax.XLabel.FontSize = fontsize;
ax.YLabel.String = 'Initial Reservoir Thickness, $d_0$ (m)';
ax.YLabel.FontSize = fontsize;
% ax.XScale = 'log';
ax.YScale = 'log';
ax.YTick = [1e-3 1e-2 1e-1 1 1e1 1e2 1e3];
ax.YTickLabel = {'$10^{-3}$','$10^{-2}$','$10^{-1}$','1','$10$','$10^{2}$','$10^{3}$'};
ax.Layer = 'Top';

cb = colorbar;
cb.FontSize = fontsize;
cb.TickLabelInterpreter = 'latex';
cb.Title.String = {'$d_{salt}$ (m)'};
cb.Title.Interpreter = 'latex';
cb.Limits = [1e-3 1e3];
cb.Ticks = [1e-3 1e-2 1e-1 1 1e1 1e2 1e3];
cb.TickLabels = {'$10^{-3}$','$10^{-2}$','$10^{-1}$','1','$10$','$10^{2}$','$10^{3}$'};

hold on
[C,h] = contour(x,y,d_salt,[1e-3 1e-2 1e-1 1 1e1 1e2 1e3],'LineColor','k');
clabel(C,h,'Interpreter','latex','FontSize',fontsize-2,'Color','k')
h.LabelSpacing = 1000;

%% Figure Formatting
%  [left bottom width height]
M = 1;
N = 2;

% left
margin = 0.175;
ax.Position(1) = margin/N*N-0.025;

% height
htot = 1-(margin/M+margin/2);
ax.Position(4) = htot/M;

% width
wtot = 1-(margin/N+2*3/4*margin);
ax.Position(3) = wtot;

% bottom
ax.Position(2) = margin/M;

f = gcf;
f.Color = 'w';
f.Units = 'inches';
width = 3.5;
height = 3;
f.Position(3:4) = [width height];

