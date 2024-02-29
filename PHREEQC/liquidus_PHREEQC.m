function [P,T,Sb] = liquidus_PHREEQC(fn)
% Extracts the temperature (C) and brine salinity (ppt) which define the
% liquidus curve for a solution composition and pressure defined by PHREEQC
% over the temperature range from the freezing point to the eutectic point.
%
% Syntax:
% [P,T,Sb] = liquidus_PHREEQC(fn)
%
% Inputs:
% fn    filename of PHREEQC output file, string
%
% Outputs:
% P     Pressure (Pa)
% T     Temperature (C)
% Sb    Brine Salinity (ppt)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%% Read PHREEQC output file
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


br_H2O = pqo.br_H2O(ind);
ms_b = pqo.br_salt(ind); % mass of salt in brine
T = pqo.T(ind); % C
mb = ms_b+br_H2O; % mass of brine
Sb = 1000*ms_b./mb;

Teut = pqo.eutectic.T;
ind = find(T>=Teut);
T = T(ind);
Sb = Sb(ind);

Tm = Tmelt(P); % K
Tm = Tm - 273.15; % C
T = [Tm; T];
Sb = [0; Sb];
end