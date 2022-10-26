function P = Pmelt(T)
% Calculates the melting pressure of ice Ih as a function of teperature
% from the ice Ih-ice III-liquid triple point to the vapor–liquid–solid
% triple point. Values outside the range of validity will output as NaN.
%
% Syntax:
% P = Pmelt(T)
%
% Inputs:
% T, Temperature    (K)
%
% Outputs:
% P, Pressure       (Pa)
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

Pt = 611.657; % Pa
Tt = 273.16; % K
a1 = 0.119539337e7;
a2 =  0.808183159e5;
a3 = 0.333826860e4;
b1 = 0.3e1;
b2 = 0.2575e2;
b3 = 0.10375e3;

theta = T/Tt;
P = Pt*(1+a1*(1-theta.^b1)+a2*(1-theta.^b2)+a3*(1-theta.^b3));

P(T > Tmax | T < Tmin) = NaN;

if any(isnan(P))
    warning(['Input temperature outside valid range (',num2str(Tmin),' K <= T <= ',num2str(Tmax),' K).'])
end
end
