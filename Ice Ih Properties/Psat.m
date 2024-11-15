function P = Psat(T)
% Calculates the saturation pressure of water as a function of teperature
% from 273.15 K to 623.15 K. Values outside the range of validity will
% output as NaN.
%
% Syntax:
% P = Psat(T)
%
% Inputs:
% T     Temperature (K), scalar or vector
%
% Outputs:
% P     Pressure (Pa), scalar or vector
%
% Source:
% IAPWS R7-97(2012)
% http://www.iapws.org/relguide/IF97-Rev.pdf
%
% Range of Validity:
% 273.15 K <= T <= 623.15 K
% 611.213 Pa <= P <= 16.521 MPa
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%%
Tmin = 273.15;
Tmax = 623.15;

Ps = 1e6; % Pa
Ts = 1; % K

n1 =  0.11670521452767e4;
n2 = -0.72421316703206e6;
n3 = -0.17073846940092e2;
n4 =  0.12020824702470e5;
n5 = -0.32325550322333e7;
n6 =  0.14915108613530e2;
n7 = -0.48232657361591e4;
n8 =  0.40511340542057e6;
n9 = -0.23855557567849; 
n10 = 0.65017534844798e3;

c = (T/Ts)+(n9./((T/Ts)-n10));

A = c.^2+n1*c+n2;
B = n3*c.^2+n4*c+n5;
C = n6*c.^2+n7*c+n8;
P = Ps*((2*C)./(-B+(B.^2-4*A.*C).^(1/2))).^4;

P(T > Tmax | T < Tmin) = NaN;

if any(isnan(P))
    warning(['Input temperature outside valid range (',num2str(Tmin),' K <= T <= ',num2str(Tmax),' K).'])
end
end