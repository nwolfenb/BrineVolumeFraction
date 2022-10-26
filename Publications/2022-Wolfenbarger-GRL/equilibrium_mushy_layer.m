% equilibrium_mushy_layer.m
% Used to generate the plots shown in Figure S2 in
%
% Wolfenbarger, N. S., Fox-Powell, M. G., Buffo, J. J., Soderlund, K. M., &
% Blankenship, D. D. (submitted). The Habitability of Brine Pockets in
% Europa's Ice Shell. Geophysical Research Letters,
% https://doi.org/10.1002/essoar.10512037.1.
%
% by following the approach of 
%
% Buffo, J. J., Schmidt, B. E., Huber, C., & Meyer, C. R. (2021).
% Characterizing the ice-ocean interface of icy worlds: A theoretical
% approach. Icarus, 360, 114318,
% https://doi.org/10.1016/j.icarus.2021.114318.
%
% and equating their equation 21 to equation 4 in
%
% Griewank, P. J., and D. Notz (2013), Insights into brine dynamics and sea
% ice desalination from a 1-D model study of gravity drainage, Journal of
% Geophysical Research: Oceans, 118, 3370–3386, doi:10.1002/jgrc.20247.

clear all; close all; clc
addpath('..\..\FREZCHEM')
addpath('..\..\Ice Ih Properties')
addpath('..\..\Water Properties')

%% Defaults
fontsize = 8;
linewidth = 2;

%% Constants
Rac = 1.01e-2;
kappa_i = 1.09e-6; % m^2/s
ci = 2e3; % J/kg/K
L = 334774; % J/kg
alpha = 1.56e-1; % kg/m^3/s
kappa_br = 1.48e-7; %m^2/s
mu = 1.88e-3; % m^2/s
phi = 0.060;
Pi = 1e-17*(1e3*phi)^3.1;

%% Europa
figure
colors = brewermap(10,'Set2');

fn = {'..\..\Freezing Simulations\FREZCHEM\v15.1\Cl_Na_Mg_SO4_seawater\Cl_Na_Mg_SO4_seawater',...
    '..\..\Freezing Simulations\FREZCHEM\v15.1\SO4_Mg_Na_Cl\SO4_Mg_Na_Cl'};
leg = {'Chloride-Dominated','Sulfate-Dominated'};

g = 1.315; % m/s^2
Ts = 100; % K
H = 10e3; % m
Soc_vec = linspace(1,100,10); % ppt

for m = 1:length(fn)
    [T_liq,Sb_liq] = liquidus_FREZCHEM(fn{m});
    [T,rho_b] = brine_density_FREZCHEM(fn{m});
    T0 = -10+273.15;
    Toc_vec = interp1(Sb_liq(Sb_liq<110),T_liq(Sb_liq<110),Soc_vec)+273.15; % K
    hint = zeros(size(Soc_vec));
    rho_boc_vec = hint;
    for n = 1:length(Soc_vec)
        Toc = Toc_vec(n);
        rho_boc = interp1(T+273.15,rho_b,Toc);
        rho_boc_vec(n) = rho_boc;
        
        N = ci*(Toc-Ts)/L;
        myfun = @(lambda,N) lambda*exp(lambda^2)*erf(lambda)-N/sqrt(pi);
        fun = @(lambda) myfun(lambda,N);
        lambda = fzero(fun,0);
        
        T_vec = linspace(T0,Toc,H*1e3);
        H0 = (Toc-T0)/(Toc-Ts)*H;
        h_vec = linspace(H0,0,H*1e3);
        rho_b_vec = interp1(T+273.15,rho_b,T_vec);
        
        delta_rho = (rho_b_vec-rho_boc)*1e3; %kg/m^3
        
        rho_34 = interp1(Sb_liq,rho_b,34);
        
        h = (((erf(lambda)*sqrt(pi)*lambda*kappa_i*exp(lambda^2))/...
            (alpha*H))+Rac)*((kappa_br*mu)./(g*delta_rho*Pi));
        
        delta_h = h-h_vec;
        ind = find(delta_h<0,1,'last' );
        hint(n) = interp1(delta_h(ind:ind+1),h(ind:ind+1),0);
        
    end
    subplot(3,1,1)
    plot(Soc_vec,hint,'o-','color',colors(m,:),...
        'markerfacecolor',colors(m,:),'linewidth',linewidth,'markersize',4)
    hold on
    
    subplot(3,1,2)
    plot(Soc_vec,Toc_vec-273.15,'o-','color',colors(m,:),...
        'markerfacecolor',colors(m,:),'linewidth',linewidth,'markersize',4)
    hold on
    
    subplot(3,1,3)
    plot(Soc_vec,rho_boc_vec-rho_34,'o-','color',colors(m,:),...
        'markerfacecolor',colors(m,:),'linewidth',linewidth,'markersize',4)
    hold on
    
end

subplot(3,1,1)
ax1 = gca;
axis tight
ax1.YLim = [0.5 2.5];
ax1.FontSize = fontsize;
ax1.TickLabelInterpreter = 'latex';
xlabel('Ocean Salinity, $S_{ocean}$ (ppt)','FontSize',fontsize,'interpreter','latex')
ylabel('Mushy Layer Thickness (m)','FontSize',fontsize,'interpreter','latex')
title(['Ice Shell Thickness: ',num2str(H/1000),' km'],'FontSize',fontsize,'interpreter','latex')


subplot(3,1,2)
ax2 = gca;
ax2.XTickLabel = [];
ax2.YDir = 'reverse';
ax2.FontSize = fontsize;
ax2.TickLabelInterpreter = 'latex';
axis tight
ylabel('Freezing Temperature, $T_f$ ($^{\circ}$C)','FontSize',fontsize,'interpreter','latex')

subplot(3,1,3)
ax3 = gca;
ax3.FontSize = fontsize;
ax3.TickLabelInterpreter = 'latex';
axis tight
xlabel('Brine Salinity, $S_{b}$ (ppt)','FontSize',fontsize,'interpreter','latex')
ylabel('Difference in Brine Density, $\Delta \rho$ (g/cm$^3$)','FontSize',fontsize,'interpreter','latex')

%% Buffo et al. (2021)
% Constants
Ts = 100; % K
Toc = 270.9; % K
Soc = 34; % ppt
Tmp = 273.15; % K
Rac = 1.01e-2;
kappa_i = 1.09e-6; % m^2/s
kappa_br = 1.48e-7; %m^2/s
ci = 2e3; % J/kg/K
L = 334774; % J/kg
alpha = 1.56e-1; % kg/m^3/s
rho_sw = 1027.347; % kg/m^3
beta = 5.836e-4; % 1/ppt
mu = 1.88e-3; % m^2/s
phi = 0.060;
Pi = 1e-17*(1e3*phi)^3.1;
Gamma = 6.6178e-2; %K kg/g
H = 10e3;
g = 1.315;

N = ci*(Toc-Ts)/L;
myfun = @(lambda,N) lambda*exp(lambda^2)*erf(lambda)-N/sqrt(pi);
fun = @(lambda) myfun(lambda,N);
lambda = fzero(fun,-1);

myfun = @(h,alpha,g,rho_sw,beta,Gamma,Tmp,Ts,H,Toc,Soc,Pi,Rac,kappa_br,lambda,kappa_i)...
    alpha*((g*rho_sw*beta*(Gamma^-1*(Tmp-(Ts+(H-h)*(Toc-Ts)/H))-Soc)*Pi*h)/...
    (kappa_br*mu)-Rac)-(erf(lambda)*sqrt(pi)*lambda*kappa_i*exp(lambda^2))/H;
fun = @(h) myfun(h,alpha,g,rho_sw,beta,Gamma,Tmp,Ts,H,Toc,Soc,Pi,Rac,kappa_br,lambda,kappa_i);
h = fzero(fun,2);

subplot(3,1,1)
plot(34,h,'ko','markerfacecolor','k','markersize',4)

subplot(3,1,2)
S = linspace(0,100);
Tf = Tmp-Gamma*S;
plot(S,Tf-273.15,'k--','markerfacecolor','k','markersize',4)
plot(34,Tmp-Gamma*34-273.15,'ko','markerfacecolor','k','markersize',4)


subplot(3,1,3)
Sint = Gamma^(-1)*(Tmp-Tf);
rho = rho_sw*beta*(Sint-34);
plot(S,rho/1000,'k--','markerfacecolor','k','markersize',4)
plot(34,0,'ko','markerfacecolor','k','markersize',4)

leg{end+1} = 'Buffo et al. (2021)';
leg{end+1} = '34 ppt';

legend(leg,'location','NorthWest','interpreter','latex','FontSize',fontsize)


%% Figure Formatting
%  [left bottom width height]
M = 3;
N = 1;

% left
margin = 0.175;
ax1.Position(1) = margin/N;
ax2.Position(1) = margin/N;
ax3.Position(1) = margin/N;

% height
htot = 1-(margin/M+margin/2);
ax1.Position(4) = htot/M;
ax2.Position(4) = htot/M;
ax3.Position(4) = htot/M;

% width
wtot = 1-(margin/N+3/4*margin);
ax1.Position(3) = wtot/N;
ax2.Position(3) = wtot/N;
ax3.Position(3) = wtot/N;

% bottom
ax3.Position(2) = margin/M-margin/M/4;
ax2.Position(2) = ax3.Position(2)+ax3.Position(4)+margin/M/4;
ax1.Position(2) = ax2.Position(2)+ax2.Position(4)+2*margin/M/2;

f = gcf;
f.Color = 'w';
f.Units = 'centimeters';
width = 95/10;
height = 3*230/4/10;
f.Position(3:4) = [width height];

% labels
annotation('textbox',[ax1.Position(1)-0.85*margin/N ax1.Position(2)+ax1.Position(4) 0.1 0.1],...
    'string','a','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')
annotation('textbox',[ax2.Position(1)-0.85*margin/N ax2.Position(2)+ax2.Position(4) 0.1 0.1],...
    'string','b','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')
annotation('textbox',[ax3.Position(1)-0.85*margin/N ax3.Position(2)+ax3.Position(4) 0.1 0.1],...
    'string','c','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')

