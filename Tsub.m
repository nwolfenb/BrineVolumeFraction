function T = Tsub(P)
% Calculates the sublimation temperature of ice Ih as a function of
% pressure from 1.9 x 10^(-40) Pa to the vapor–liquid–solid triple point.
% Values outside the valid range will output as NaN. 
%
% Syntax:
% T = Tsub(P)
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
% 50 <= T <= 273.16
% 1.9 x 10^(-40) Pa <= P <= 611.6570 Pa
%
% Author:
% Natalie Wolfenbarger
% nwolfenb@utexas.edu
%
%%
Tmin = 50;
Tmax = 273.16;
Tmod = linspace(Tmin,Tmax,1e6);
Pmin = Psub(Tmin);
Pmax = Psub(Tmax);

Pt = 611.657; % Pa
Tt = 273.16; % K
a1 = -0.212144006e2;
a2 =  0.273203819e2;
a3 = -0.610598130e1;
b1 = 0.333333333e-2;
b2 = 0.120666667e1;
b3 = 0.170333333e1;

theta = Tmod/Tt;
Pmod = Pt*exp((1./theta).*(a1*theta.^b1+a2*theta.^b2+a3*theta.^b3));

T = interp1(Pmod,Tmod,P);

T(P < Pmin | P > Pmax) = NaN;

if any(isnan(T))
    warning(['Input pressure outside valid range (',num2str(Pmin,'%1.3e'),...
        ' Pa <= P <= ',num2str(Pmax,'%3.3f'),' Pa).'])
end
end