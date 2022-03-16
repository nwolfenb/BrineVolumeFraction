function T = Tmelt(p)
% Calculates the melting temperature for Ice Ih given pressure (Pa) from
% the ice Ih-lice III-liquid triple point to the vapor–liquid–solid triple
% point. Pressures outside the valid range will output as NaN.
%
% Syntax:
% T = Tmelt(p)
%
% Inputs:
% Pressure      (Pa)
%
% Outputs:
% Temperature   (K)
%
% Source:
% https://doi.org/10.1063/1.3657937
%
% Range of Validity:
% 611.6570 Pa <= P <= 208.5666 MPa
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%%
    Tmod = linspace(251.165,273.16,1e6);
    
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
    
    T = interp1(Pmod,Tmod,p); 
    
    T(p > 208.5666e6 | p < 611.6570) = NaN;
    
if ~all(p <= 208.5666e6 & p >= 611.6570)
    warning('Input pressure outside valid range (611.6570 Pa <= P <= 208.5666 MPa).')
end
end