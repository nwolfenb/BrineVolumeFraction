function [P,T,Sb] = liquidus_FREZCHEM(fn)
% Extracts the temperature (C) and brine salinity (ppt) which define the
% liquidus curve for a solution composition and pressure defined by
% FREZCHEM version 13.3 or 15.1 over the temperature range from the
% freezing point to the eutectic point.
%
% Syntax:
% [P,T,Sb] = liquidus_FREZCHEM(fn)
%
% Inputs:
% fn    filename of FREZCHEM output file, string
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
%% Read FREZCHEM output file
FrOut = read_FrOut(fn);

%% Check Pressure
p = FrOut.pressure; % bar
if ~all(p==p(1))
    warning('Pressure appears to change between temperature steps in FREZCHEM output file, adopting value at initial temperature step.')
end
P=p(1)*100000; % Pa;

pt = 611.657; % Pa, Triple point
if P < pt
    error('Input pressure must be above 611.657 Pa for brine to be stable.')
end

%% Extract Liquidus
ice_H2O = FrOut.ice_H2O;
ind = find(ice_H2O>0); % eliminate output above freezing point

br_H2O = FrOut.br_H2O(ind);
ms_b = FrOut.br_salt(ind); % mass of salt in brine
T = FrOut.T(ind) - 273.15; % K to C
mb = ms_b+br_H2O; % mass of brine
Sb = 1000*ms_b./mb;

Tm = Tmelt(P); % K
Tm = Tm - 273.15; % C
T = [Tm; T];
Sb = [0; Sb];
end