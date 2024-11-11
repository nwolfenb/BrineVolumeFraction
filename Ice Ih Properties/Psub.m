function P = Psub(T)
% Calculates the sublimation pressure of ice Ih as a function of teperature
% from the vapor–liquid–solid triple point to 50 K. Values outside the
% range of validity will output as NaN.
%
% Syntax:
% P = Psub(T)
%
% Inputs:
% T     Temperature (K), scalar or vector
%
% Outputs:
% P     Pressure (Pa), scalar or vector
%
% Source:
% https://doi.org/10.1063/1.3657937
% IAPWS R14-08(2011)
%
% Range of Validity:
% 50 K <= T <= 273.16 K
% 1.9 x 10^(-40) Pa <= P <= 611.6570 Pa
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%%
Tmin = 50;
Tmax = 273.16;

Pt = 611.657; % Pa
Tt = 273.16; % K
a1 = -0.212144006e2;
a2 =  0.273203819e2;
a3 = -0.610598130e1;
b1 = 0.333333333e-2;
b2 = 0.120666667e1;
b3 = 0.170333333e1;

theta = T/Tt;
P = Pt*exp((1./theta).*(a1*theta.^b1+a2*theta.^b2+a3*theta.^b3));

P(T > Tmax | T < Tmin) = NaN;

if any(isnan(P))
    warning(['Input temperature outside valid range (',num2str(Tmin),' K <= T <= ',num2str(Tmax),' K).'])
end
end