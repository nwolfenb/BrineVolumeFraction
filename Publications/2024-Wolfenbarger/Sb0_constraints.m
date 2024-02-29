% Sb0_constraints.m 
% Used to generate the plot shown in Figure 4 in
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
rhob = interp1(Sb,rho_b,S);
d0 = logspace(-1,4,1e3)';
d_salt = d0*S/1000.*(rhob/rhos)*mm_MgSO411H2O/mm_MgSO4*(1/veut);

X = d_salt;
Y = repmat(d0,1,length(S));
V = repmat(S,length(d0),1);

x = X(:);
y = Y(:);
v = V(:);

xq = logspace(-1,4,1e3);
yq = logspace(-1,4,1e3);

Xq = repmat(xq,length(yq),1);
Yq = repmat(yq',1,length(xq));
F = scatteredInterpolant(x,y,v,'linear','none');
vq = F(Xq(:),Yq(:));
Vq = reshape(vq',length(yq),length(xq));

figure
% Extra subplots
subplot(2,2,1)
semilogx([7.9 7.9],[0 1],'k:','LineWidth',linewidth/2) % VHF
hold on
semilogx([79.1 79.1],[0 1],'k:','LineWidth',linewidth/2) % HF
xlim([min(xq) max(xq)])
ylim([0 1])
axis off
ax11 = gca;
ax11.XScale = 'log';


subplot(2,1,2)
p = pcolor(Xq,Yq,Vq);
p.LineStyle = 'none';
p.FaceColor = 'interp';
clim([1 Seut])
axis tight

ax1 = gca;
ax1.ColorScale = 'log';
ax1.FontSize = fontsize;
ax1.XLabel.String = 'Salt Layer Thickness, $d_{salt}$ (m)';
ax1.XLabel.FontSize = fontsize;
ax1.YLabel.String = 'Initial Reservoir Thickness, $d_0$ (m)';
ax1.YLabel.FontSize = fontsize;
ax1.XScale = 'log';
ax1.YScale = 'log';
ax1.XTick = [1e-1 1 1e1 1e2 1e3 1e4];
ax1.YTick = [1e-1 1 1e1 1e2 1e3 1e4];
ax1.XTickLabel = {'$10^{-1}$','1','$10$','$10^{2}$','$10^{3}$','$10^{4}$'};
ax1.YTickLabel = {'$10^{-1}$','1','$10$','$10^{2}$','$10^{3}$','$10^{4}$'};
ax1.Layer = 'Top';

cb = colorbar;
cb.FontSize = fontsize;
cb.TickLabelInterpreter = 'latex';
cb.Title.String = {'$S_{b,0}$ (ppt)'};
cb.Title.Interpreter = 'latex';
cb.Limits = [1 Seut];
cb.Ticks = [1 10 100];
cb.TickLabels = {'1','10','100'};

keq = 0.07;
hold on
[C,h] = contour(Xq,Yq,Vq,round(Seut*keq)*ones(1,2),'LineColor','k');
clabel(C,h,'Interpreter','latex','FontSize',fontsize-2,'Color','k')


%% Below Range Resolution
patch([0.7 7.9 7.9 0.7],[0.1 0.1 1e4 1e4],100/255*ones(1,3),'EdgeColor','none')
text(exp((log(0.7)+log(7.9))/2),exp((log(0.1)+log(1e4))/2),'Unresolvable',...
    'FontSize',fontsize,'Rotation',90,'Interpreter','tex',...
    'HorizontalAlignment','Center','VerticalAlignment','Middle',...
    'FontWeight','Normal','Color','w')

patch([0.1 0.7 0.7 0.1],[0.1 0.1 1e4 1e4],0*ones(1,3),'EdgeColor','none')
text(exp((log(0.1)+log(0.7))/2),exp((log(0.1)+log(1e4))/2),'Undetectable',...
    'FontSize',fontsize,'Rotation',90,'Interpreter','tex',...
    'HorizontalAlignment','Center','VerticalAlignment','Middle',...
    'FontWeight','Normal','Color','w')

plot([7.9 7.9],[0.1 1e4],'k:','LineWidth',linewidth/2)
text(7.9,1e4,{'REASON     ';''},...
    'FontSize',fontsize-2,'Interpreter','tex',...
    'HorizontalAlignment','right','VerticalAlignment','Bottom')
text(7.9,1e4,{'VHF';''},...
    'FontSize',fontsize-2,'Interpreter','tex',...
    'HorizontalAlignment','center','VerticalAlignment','Bottom')

plot([28.2 28.2],[0.1 1e4],'k:','LineWidth',linewidth/2)
text(28.2,1e4,'RIME',...
    'FontSize',fontsize-2,'Interpreter','tex',...
    'HorizontalAlignment','center','VerticalAlignment','Bottom')

plot([79.1 79.1],[0.1 1e4],'k:','LineWidth',linewidth/2)
text(79.1,1e4,{'HF';''},...
    'FontSize',fontsize-2,'Interpreter','tex',...
    'HorizontalAlignment','center','VerticalAlignment','Bottom')

text(exp((log(7.9)+log(1e4))/2),exp((log(8.5)+log(1e4))/2),'Ocean Injection',...
    'FontSize',fontsize,'Rotation',45,'Interpreter','tex',...
    'HorizontalAlignment','Center','VerticalAlignment','Bottom',...
    'FontWeight','Bold','Color','w')


%% Figure Formatting
%  [left bottom width height]
M = 1;
N = 2;

% left
margin = 0.175;
ax1.Position(1) = margin/N*N-0.025;
ax11.Position(1) = ax1.Position(1);

% height
htot = 1-(margin/M+margin/2);
ax1.Position(4) = htot/M;
ax11.Position(4) = margin/N/2;

% width
wtot = 1-(margin/N+2*3/4*margin);
ax1.Position(3) = wtot;
ax11.Position(3) = ax1.Position(3);

% bottom
ax1.Position(2) = margin/M;
ax11.Position(2) = ax1.Position(2)+ax1.Position(4);

f = gcf;
f.Color = 'w';
f.Units = 'inches';
width = 3.5;
height = 3;
f.Position(3:4) = [width height];
