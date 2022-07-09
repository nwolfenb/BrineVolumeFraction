function T = Tsat(P)
% Calculates the saturation pressure of water as a function of pressure
% from 611.213 Pa to 16.521 MPa. Values outside the valid range will output
% as NaN.
%
% Syntax:
% T = Tsat(P)
%
% Inputs:
% P, Pressure       (Pa)
%
% Outputs:
% T, Temperature    (K)
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
% nwolfenb@utexas.edu
%
%%
Tmin = 273.15;
Tmax = 623.15;
Tmod = linspace(Tmin,Tmax,1e6);
Pmin = Psat(Tmin);
Pmax = Psat(Tmax);

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

c = (Tmod/Ts)+(n9./((Tmod/Ts)-n10));

A = c.^2+n1*c+n2;
B = n3*c.^2+n4*c+n5;
C = n6*c.^2+n7*c+n8;
Pmod = Ps*((2*C)./(-B+(B.^2-4*A.*C).^(1/2))).^4;

T = interp1(Pmod,Tmod,P);

T(P < Pmin | P > Pmax) = NaN;

if any(isnan(T))
    warning(['Input pressure outside valid range (',num2str(Pmin,'%1.3e'),...
        ' Pa <= P <= ',num2str(Pmax,'%3.3f'),' Pa).'])
end
end