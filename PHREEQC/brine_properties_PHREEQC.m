function [P,T,Sb,rho_b,sigma] = brine_properties_PHREEQC(fn)
% Extracts the temperature (C), salinity (ppt) and brine density (g/cm^3)
% for a solution composition and pressure defined by PHREEQC over the
% temperature range from the freezing point to the eutectic point.
%
% Syntax:
% [P,T,Sb,rho_b,sigma] = brine_properties_PHREEQC(fn)
%
% Inputs:
% fn    Filename of PHREEQC output file, string
%
% Outputs:
% P      Pressure (Pa), scalar
% T      Temperature (C), vector
% Sb     Brine Salinity (ppt), vector
% rho_b  Brine Density (g/cm^3), vector
% sigma  Brine Conductivity (uS/cm), vector
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%% Read PHREEQC output file
pqo = read_pqo(fn);

%% Check Pressure
p = pqo.pressure; % atm
if ~all(p==p(1))
    warning('Pressure appears to change between temperature steps in FREZCHEM output file, adopting value at initial temperature step.')
end
P=p(1)*101325; % Pa;

pt = 611.657; % Pa, Triple point
if P < pt
    error('Input pressure must be above 611.657 Pa for brine to be stable.')
end

%% Extract Brine Properties
ice_H2O = pqo.ice_H2O;
ind = find(ice_H2O>0 & pqo.T>=pqo.eutectic.T); % eliminate output above freezing point and below eutectic

br_H2O = pqo.br_H2O(ind);
ms_b = pqo.br_salt(ind); % mass of salt in brine
T = pqo.T(ind); % C
mb = ms_b+br_H2O; % mass of brine
Sb = 1000*ms_b./mb;
rho_b = pqo.rho_b(ind);
sigma = pqo.cond(ind);


Tm = Tmelt(P); % K
Tm = Tm - 273.15; % C
T = [Tm; T];
Sb = [0; Sb];
rho_b = [water_density(Tm+273.15,P)/1e3; rho_b];
sigma = [0; sigma];
end