function T = Tmelt(P)
% Calculates the melting temperature of ice Ih as a function of pressure
% from the ice Ih-ice III-liquid triple point to the vapor–liquid–solid
% triple point. Values outside the valid range will output as NaN. 
%
% Syntax:
% T = Tmelt(p)
%
% Inputs:
% P, Pressure       (Pa)
%
% Outputs:
% T, Temperature    (K)
%
% Source:
% https://doi.org/10.1063/1.3657937
% IAPWS R14-08(2011)
%
% Range of Validity:
% 251.165 K <= T <= 273.16 K
% 611.657 Pa <= P <= 208.567 MPa
%
% Author:
% Natalie Wolfenbarger
% nwolfenb@utexas.edu
%
%%
Tmin = 251.165;
Tmax = 273.16;
Tmod = linspace(Tmin,Tmax,1e6);
Pmin = Pmelt(Tmax);
Pmax = Pmelt(Tmin);

Pt = 611.657; % Pa
Tt = 273.16; % K
a1 = 0.119539337e7;
a2 =  0.808183159e5;
a3 = 0.333826860e4;
b1 = 0.3e1;
b2 = 0.2575e2;
b3 = 0.10375e3;

theta = Tmod/Tt;
Pmod = Pt*(1+a1*(1-theta.^b1)+a2*(1-theta.^b2)+a3*(1-theta.^b3));

T = interp1(Pmod,Tmod,P);

T(P < Pmin | P > Pmax) = NaN;

if any(isnan(T))
    warning(['Input pressure outside valid range (',num2str(Pmin,'%3.3f'),...
        ' Pa <= P <= ',num2str(Pmax/1e6,'%3.3f'),' MPa).'])
end
end