function [P,T,rho_b] = brine_density_PHREEQC(fn)
% Extracts the pressure (Pa), temperature (C) and brine density (g/cm^3)
% for a solution composition and pressure defined by the PHREEQC output
% over the temperature range from the freezing point to the eutectic point.
%
% Inputs:
% fn    filename of PHREEQC output file, string
%
% Outputs:
% P     Pressure (MPa)
% T     Temperature (C)
% rho_b  Brine Density (g/cm^3)

pqo = read_pqo(fn);

%% Check Pressure
p = pqo.pressure; % bar
if ~all(p==p(1))
    if all(isnan(p))
        P = 101325; % Pa
    else
        warning('Pressure appears to change between temperature steps in PHREEQC output file, adopting value at initial temperature step.')
    end
else
P=p(1)*101325; % Pa
end

pt = 611.657; % Pa, Triple point
if P < pt
    error('Input pressure must be above 611.657 Pa for brine to be stable.')
end

%% Extract Liquidus
ice_H2O = pqo.ice_H2O;
ind = find(ice_H2O>0); % eliminate output above freezing point

T = pqo.T(ind); % C
rho_b = pqo.rho_b(ind);

Teut = pqo.eutectic.T;
ind = find(T>=Teut);
T = T(ind);
rho_b = rho_b(ind);

Tm = Tmelt(P); % K
Tm = Tm - 273.15; % C
T = [Tm; T];
rho_b = [water_density(Tm+273.15,P)/1e3; rho_b];
end