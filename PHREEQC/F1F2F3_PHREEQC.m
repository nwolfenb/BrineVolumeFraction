function [P,T,F1,F2,F3,k_star,rho_ss] = F1F2F3_PHREEQC(fn)
% Calculates the phase behavior functions F1 and F2 in Cox and Weeks (1983)
% and a new phase behavior function F3 from freezing simulations in PHREEQC
% over a temperature range from the pure ice pressure melting point to the
% eutectic point.
%
% Inputs:
% fn    filename of FREZCHEM output file, string
%
% Outputs:
% P         pressure specified in freezing simulation
% T         temperature range (pure ice pressure melting point to eutectic)
% F1        rho_b*Sb*(1+k), k = ms_ss/ms_b
% F2        (1+C)*rho_b/rho_i-C*rho_b/rho_ss-1), C = mss/mb 
% F3        C*rho_b/rho_ss, C = mss/mb
% k_star    scale factor to account for hydration of minerals at/beyond eutectic
% rho_ss    total density of solid salt precipitating at/beyond the eutectic

%% Check if melt is stable for simulated pressure
pqo = read_pqo(fn);
p = pqo.pressure; % atm
if ~all(p==p(1))
    warning('Pressure appears to change between temperature steps in FREZCHEM output file, adopting value at initial temperature step.')
end
P=p(1)*101325; % Pa;

pt = 611.657; % Pa, Triple point
if P < pt
    error('Brine is not stable for the pressure specified in the FREZCHEM file, pressure must be greater than the triple point (611.657 Pa = 0.00611657 bar)')
else
    
    %% Extract Data from PHREEQC file
    ice_H2O = pqo.ice_H2O;
    ms_b0 = pqo.br_salt(1); % initial mass of salts in brine
    
    ind = find(ice_H2O>0 & pqo.T>pqo.eutectic.T); % eliminate output above freezing point and below eutectic
    
    br_H2O = pqo.br_H2O(ind);
    ms_b = pqo.br_salt(ind); % mass of salt in brine
    mss = pqo.sm(ind); % mass of solid salts
    rho_b = pqo.rho_b(ind);
    rho_ss = pqo.rho_ss(ind);
    
    T_pqo = pqo.T(ind); % C
    
    %% Defining Parameters for Brine Volume Calculations
    mb = ms_b+br_H2O; % mass of brine
    ms_ss = ms_b0-ms_b; % mass of salts not in brine
    k = ms_ss./ms_b;
    C = mss./mb;
    Sb = 1000*ms_b./mb;
    rho_i = ice_density(T_pqo+273.15,P)/1e3;
    
    %% F1, F2, and F3
    F1 = rho_b.*Sb.*(1+k);
    F2 = (1+C).*(rho_b./rho_i)-(C.*rho_b./rho_ss)-1;
    F3 = C.*rho_b./rho_ss;
    ind = find(isnan(rho_ss));
    F2(ind) = rho_b(ind)./rho_i(ind)-1;
    F3(ind) = 0;
    
    %% Salts forming at/beyond the eutectic
    k_star = max(pqo.sm)/ms_b0;
    rho_ss = max(pqo.rho_ss);

    
    %% Melting Point
    Tm = Tmelt(P); % K
    Tm = Tm - 273.15; % C
    
    %% Extension of F1, F2, and F3 to T = Tm
    T = [Tm; T_pqo];
    F1 = [0; F1];
    F2 = [water_density(Tm+273.15,P)/ice_density(Tm+273.15,P)-1; F2];
    F3 = [0; F3];
    
end    
end