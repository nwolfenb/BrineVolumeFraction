function [P,T,Sb,rho_b] = brine_properties_FREZCHEM(fn)
% Extracts the temperature (C), salinity (ppt) and brine density (g/cm^3)
% for a solution composition and pressure defined by FREZCHEM version 13.3
% or 15.1 over the temperature range from the freezing point to the
% eutectic point.
%
% Syntax:
% [P,T,Sb,rho_b] = brine_properties_FREZCHEM(fn)
%
% Inputs:
% fn    Filename of FREZCHEM output file, string
%
% Outputs:
% P     Pressure (Pa), scalar
% T     Temperature (C), vector
% Sb    Brine Salinity (ppt), vector
% rho_b Brine Density (g/cm^3), vector
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%% Read FREZCHEM output file
FrOut = read_FrOut(fn);

%% Check Pressure
p = FrOut.pressure; % bar
if ~all(p==p(1))
    warning('Pressure appears to change between temperature steps in FREZCHEM output file, adopting value at initial temperature step.')
end
P=p(1)*100000; % Pa;

Pt = 611.657; % Pa, Triple point
if P < Pt
    error('Input pressure must be above 611.657 Pa for brine to be stable.')
end

%% Extract Brine Properties
ice_H2O = FrOut.ice_H2O;
ind = find(ice_H2O>0); % eliminate output above freezing point

br_H2O = FrOut.br_H2O(ind);
ms_b = FrOut.br_salt(ind); % mass of salt in brine
T = FrOut.T(ind) - 273.15; % K to C
mb = ms_b+br_H2O; % mass of brine
Sb = 1000*ms_b./mb;
rho_b = FrOut.rho_b(ind);

Tm = Tmelt(P); % K
Tm = Tm - 273.15; % C
T = [Tm; T];
Sb = [0; Sb];
rho_b = [water_density(Tm+273.15,P)/1e3; rho_b];
end