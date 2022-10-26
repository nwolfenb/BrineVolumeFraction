function [T,rhob] = brine_density_FREZCHEM(fn)
% Extracts the temperature (C) and brine density (g/cm^3) for a solution
% composition and pressure defined by FREZCHEM version 13.3 or 15.1 over
% the temperature range from the freezing point to the eutectic point.
%
% Syntax:
% [T,rhob] = brine_density_FREZCHEM(fn)
%
% Inputs:
% fn    filename of FREZCHEM output file, string
%
% Outputs:
% T     Temperature (C)
% rhob  Brine Density (g/cm^3)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%

%% Read FREZCHEM file
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

%% Extract Liquidus
ice_H2O = FrOut.ice_H2O;
ind = find(ice_H2O>0); % eliminate output above freezing point

T = FrOut.T(ind)-273.15; % C
rhob = FrOut.rho_b(ind);

Tm = Tmelt(P); % K
Tm = Tm - 273.15; % C
T = [Tm; T];
rhob = [water_density(Tm+273.15,P)/1e3; rhob];
end