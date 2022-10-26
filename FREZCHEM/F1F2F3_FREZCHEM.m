function [T,F1,F2,F3] = F1F2F3_FREZCHEM(fn)
% Calculates the phase behavior functions F1 and F2 in Cox and Weeks (1983)
% and a new phase behavior function F3 for a solution composition and
% pressure defined by FREZCHEM version 13.3 or 15.1 over the temperature
% range from the pure ice pressure melting point to the eutectic point.
%
% Inputs:
% fn    filename of FREZCHEM output file, string
%
% Outputs:
% T     temperature range (pure ice pressure melting point to eutectic)
% F1    rho_b*Sb*(1+k), k = ms_ss/ms_b
% F2    (1+C)*rho_b/rho_i-C*rho_b/rho_ss-1), C = mss/mb 
% F3     C*rho_b/rho_ss, C = mss/mb


%% Check if melt is stable for simulated pressure
FrOut = read_FrOut(fn);
p = FrOut.pressure; % bar
if ~all(p==p(1))
    warning('Pressure appears to change between temperature steps in FREZCHEM output file, adopting value at initial temperature step.') 
end
p=p(1)*100000; % Pa;

pt = 611.657; % Pa, Triple point
if p < pt
    error('Brine is not stable for the pressure specified in the FREZCHEM file, pressure must be greater than the triple point (611.657 Pa = 0.00611657 bar)')
else
    
    %% Extract Data from FREZCHEM file
    ice_H2O = FrOut.ice_H2O;
    S0 = FrOut.br_salt(1); % initial mass of salts in brine
    
    ind = find(ice_H2O>0); % eliminate output above freezing point
    
    br_H2O = FrOut.br_H2O(ind);
    ms_b = FrOut.br_salt(ind); % mass of salt in brine
    mss = FrOut.sm(ind); % mass of solid salts
    rho_b = FrOut.rho_b(ind);
    rho_ss = FrOut.rho_ss(ind);
    
    T_FrOut = FrOut.T(ind) - 273.15; % K to C
    
    %% Defining Parameters for Brine Volume Calculations
    mb = ms_b+br_H2O; % mass of brine
    ms_ss = S0-ms_b; % mass of salts not in brine
    k = ms_ss./ms_b;
    C = mss./mb;
    Sb = 1000*ms_b./mb;
    rho_i = ice_density(T_FrOut+273.15,p)/1e3;
    
    %% F1, F2, and F3
    F1 = rho_b.*Sb.*(1+k);
    F2 = (1+C).*(rho_b./rho_i)-(C.*rho_b./rho_ss)-1;
    F3 = C.*rho_b./rho_ss;
    ind = find(isnan(rho_ss));
    F2(ind) = rho_b(ind)./rho_i(ind)-1;
    F3(ind) = 0;
    
    %% Melting Point
    Tm = Tmelt(p); % K
    Tm = Tm - 273.15; % C
    
    %% Extension of F1, F2, and F3 to T = Tm
    T = [Tm; T_FrOut];
    F1 = [0; F1];
    F2 = [water_density(Tm+273.15,p)/ice_density(Tm+273.15,p)-1; F2];
    F3 = [0; F3];
    
end    

end