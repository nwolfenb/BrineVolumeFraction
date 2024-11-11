function [P,Tmat,Smat,Vi_V, Vb_V, Vs_V] = volume_fraction_FREZCHEM(T,S,fn)
% Calculates volume fractions of each component (ice, brine, solid salts)
% for an impure ice with salinity, S, temperature, T, for a solution
% composition and pressure defined by FREZCHEM version 13.3 or 15.1 over
% the temperature range from the freezing point to the eutectic point. If T
% and S are both vectors, the brine volume fraction is calculated for each
% salinity and temperature where temperature varies across columns and
% salinity varies across rows.
%
% Syntax:
% [P,Tmat,Smat,Vi_V, Vb_V, Vs_V] = volume_fraction_FREZCHEM(T,S,fn)
%
% Inputs:
% T     Temperature (C), scalar or vector
% S     Bulk salinity of ice (ppt), scalar or vector
% fn    Filename of FREZCHEM output file, string
%
% Outputs:
% P     Pressure (Pa), scalar
% Tmat  Temperature (C), vector or matrix
% Smat  Bulk salinity of ice (ppt), vector or matrix
% Vi_V  Ice volume fraction, vector or matrix
% Vb_V  Brine volume fraction, vector or matrix
% Vs_V  Salt volume fraction, vector or matrix
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%% Check Inputs
if ~iscolumn(T)
    T = T'; % column vector
end

if iscolumn(S)
    S = S'; % row vector
end

%% Check if melt is stable for simulated pressure
FrOut = read_FrOut(fn);
p = FrOut.pressure; % bar
if ~all(p==p(1))
    warning('Pressure appears to change between temperature steps in FREZCHEM output file, adopting value at initial temperature step.')
end
P=p(1)*100000; % Pa;

pt = 611.657; % Pa, Triple point
if P < pt
    Vb_V = zeros(size(T));
    
    fw = S/1000; % assumes anhydrous salt mass = hydrated salt mass
    rho_ss = 1.5; % average solid salt density (Cox and Weeks, 1983)
    rho_i = ice_density(T,P);
    Vs_V = fw/(fw+(1-fw)*(rho_ss/rho_i));
    
    Vi_V = 1-Vs_V;
    
    Vb_V(T>(Tsub(P)-273.15)) = NaN;
    Vi_V(T>(Tsub(P)-273.15)) = NaN;
    Vs_V(T>(Tsub(P)-273.15)) = NaN;
    warning(['Input pressure must be above 611.657 Pa for brine to be stable.',...
        ' Volume fractions of ice and salt estimated using the salinity',...
        ' of the ice assuming a salt density of 1.5 g/cm^3.'])
    % Use anhydrous salt mass to estimate average anhydrous salt density to
    % estimate volume fractions
else
    
    %% Extract Data from FREZCHEM file
    ice_H2O = FrOut.ice_H2O;
    ms_b0 = FrOut.br_salt(1); % initial mass of salts in brine
    
    ind = find(ice_H2O>0); % eliminate output above freezing point
    
    br_H2O = FrOut.br_H2O(ind);
    ms_b = FrOut.br_salt(ind); % mass of salt in brine
    mss = FrOut.sm(ind); % mass of solid salts
    rho_b = FrOut.rho_b(ind);
    rho_ss = FrOut.rho_ss(ind);
    
    T_FrOut = FrOut.T(ind) - 273.15; % K to C
    
    %% Defining Parameters for Brine Volume Calculations
    mb = ms_b+br_H2O; % mass of brine
    ms_ss = ms_b0-ms_b; % mass of salts not in brine
    k = ms_ss./ms_b;
    C = mss./mb;
    Sb = 1000*ms_b./mb;
    rho_i = ice_density(T_FrOut+273.15,P)/1e3;
    
    %% F1 and F2
    F1 = rho_b.*Sb.*(1+k);
    F2 = (1+C).*(rho_b./rho_i)-(C.*rho_b./rho_ss)-1;
    ind = find(isnan(rho_ss));
    F2(ind) = rho_b(ind)./rho_i(ind)-1;
    
    %% Freezing Point
    Tm = Tmelt(P); % K
    Tm = Tm - 273.15; % C
    
    %% Salt Volume
    % You know where Vb_V>1 and Vb_V<0, ice is not forming; however, this
    % does not mean salt is not forming. Where this occurs, you need to
    % solve the system:
    % Vs_V = C.*(rho_b./rho_ss).*Vb_V,
    % Vs_V + Vb_V = 1 
    
    Vb_V = rho_i*S./(F1-rho_i.*F2*S);
    F3 = C.*(rho_b./rho_ss);
    F3 = repmat(F3,1,length(S));
    Vs_V = F3./(1+F3);
    
    ind_bad = find(Vb_V+Vs_V>1 | Vb_V<0 | Vb_V>1);
    ind_good = find(Vb_V+Vs_V<=1 & Vb_V>=0);
    Vb_V(ind_bad) = 1 - Vs_V(ind_bad);
    Vs_V(ind_good) = F3(ind_good).*Vb_V(ind_good);
    
    T_FrOut_mod = [Tm; T_FrOut];
    Vs_V_max = max(Vs_V);
    Vs_V_max = repmat(Vs_V_max,length(T),1);
    Vs_V = [zeros(1,length(S)); Vs_V];
    Tmat = repmat(T,1,length(S));
    Smat = repmat(S,length(T),1);
    if isscalar(S)
        Vs_V = interp1(T_FrOut_mod,Vs_V,T,'linear');
    else
        [Smat_FrOut,Tmat_FrOut] = meshgrid(S,T_FrOut_mod);
        Vs_V = interp2(Smat_FrOut,Tmat_FrOut,Vs_V,Smat,Tmat,'linear');
        if C(1)>0
            Vs_V(Tmat(:,1)>T_FrOut(1),:) = 0;
        end
    end
    
    if isnan(FrOut.eutectic.sm)
        Vs_V(Tmat<T_FrOut(end)) = Vs_V_max(Tmat<T_FrOut(end));
        warning(['Salt volume fraction below the eutectic temeprature ',...
            'is set equal to the salt volume fraction at the temperature',...
            ' step immediately prior to the eutctic.'])
        if all(Vs_V==0,'all')
            Vs_V(Tmat<T_FrOut(end)) = NaN;
            warning(['Salt volume fraction below the euectic could not be',...
                ' defined. If mineral precipitating at the eutectic',...
                ' is known, Vs_V can be calculated manually.'])
        end
    else
        rho_ss = FrOut.eutectic.rho_ss*1000;
        rho_i = ice_density(T(T<T_FrOut(end))+273.15,P);
        if ~isempty(rho_i)
            rho_i = repmat(rho_i,1,length(S));
            
            k_star = FrOut.eutectic.sm/ms_b0;
            Vs_V(Tmat<T_FrOut(end)) = S*k_star./...
                (S*k_star+(1000-S*k_star).*(rho_ss./rho_i));
            
            disp(['k* = ',num2str(k_star)])
        end
    end
    
    %% Extension of F1 and F2 to T = Tm
    T_FrOut = [Tm; T_FrOut];
    F1 = [0; F1];
    F2 = [water_density(Tm+273.15,P)/ice_density(Tm+273.15,P)-1; F2];
    
    F1 = interp1(T_FrOut,F1,T,'linear');
    F2 = interp1(T_FrOut,F2,T,'linear');
    rho_i = ice_density(T+273.15,P)/1e3;
    
    %% Brine Volume
    Vb_V = rho_i*S./(F1-rho_i.*F2*S);

    ind = find(Vb_V+Vs_V>1 | Vb_V<0 | Vb_V>1);
    Vb_V(ind) = 1 - Vs_V(ind);
    Vb_V(isnan(Vb_V)) = 0;
    Vs_V(isnan(Vs_V)) = 0;
    Vb_V(Vb_V==0 & Vs_V==0 & Smat>0) = 1;
    
    %% Ice Volume
    Vi_V = 1-(Vb_V+Vs_V);
    
end
end