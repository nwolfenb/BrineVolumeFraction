function rho = water_density(T,P)
% Density of water as a function of temperature and pressure from the Gibbs
% energy equation of state.
%
% Syntax:
% rho = water_density(T,P)
%
% Inputs:
% T, Temperature    (K)
% P, Pressure       (Pa)
%
% Outputs:
% rho, Density      (kg/m^3)
%
% Source:
% IAPWS SR7-09(2009)
% http://www.iapws.org/relguide/OceanLiquid.pdf
%
% Range of Validity:
% 100 Pa <= P <= 100 MPa
% (270.5 – P × 7.43 × 10–8 Pa–1) K <= T <= 313.15 K
%
% Author:
% Natalie Wolfenbarger
% nwolfenb@utexas.edu
%
%% Convert to Kelvin
if min(T)<=0
    T = T+273.15;
end

%% Coefficients of the equation of state (Gibbs potential function)
g = zeros(8,7);
g(1,1) = 0.101342743139674e3;
g(1,2) = 0.100015695367145e6;
g(1,3) = -0.254457654203630e4;
g(1,4) = 0.284517778446287e3;
g(1,5) = -0.333146754253611e2;
g(1,6) = 0.420263108803084e1;
g(1,7) = -0.546428511471039;
g(2,1) = 0.590578347909402e1;
g(2,2) = -0.270983805184062e3;
g(2,3) = 0.776153611613101e3;
g(2,4) = -0.196512550881220e3;
g(2,5) = 0.289796526294175e2;
g(2,6) = -0.213290083518327e1;
g(3,1) = -0.123577859330390e5;
g(3,2) = 0.145503645404680e4;
g(3,3) = -0.756558385769359e3;
g(3,4) = 0.273479662323528e3;
g(3,5) = -0.555604063817218e2;
g(3,6) = 0.434420671917197e1;
g(4,1) = 0.736741204151612e3;
g(4,2) = -0.672507783145070e3;
g(4,3) = 0.499360390819152e3;
g(4,4) = -0.239545330654412e3;
g(4,5) = 0.488012518593872e2;
g(4,6) = -0.166307106208905e1;
g(5,1) = -0.148185936433658e3;
g(5,2) = 0.397968445406972e3;
g(5,3) = -0.301815380621876e3;
g(5,4) = 0.152196371733841e3;
g(5,5) = -0.263748377232802e2;
g(6,1) = 0.580259125842571e2;
g(6,2) = -0.194618310617595e3;
g(6,3) = 0.120520654902025e3;
g(6,4) = -0.552723052340152e2;
g(6,5) = 0.648190668077221e1;
g(7,1) = -0.189843846514172e2;
g(7,2) = 0.635113936641785e2;
g(7,3) = -0.222897317140459e2;
g(7,4) = 0.817060541818112e1;
g(8,1) = 0.305081646487967e1;
g(8,2) = -0.963108119393062e1;

%% Constants
T0 = 273.15; % K
Ts = 40; % K
p0 = 101325; % K
ps = 1e8; % K
gs = 1; % J/kg

%% Dimensionless variables
tau = (T-T0)/Ts;
pii = (P-p0)/ps;

%% gp
gsum = 0;
for j = 0:7
    for k = 1:6
        k_ind = k+1;
        j_ind = j+1;
        gsum = gsum + (k*g(j_ind,k_ind)*tau.^(j).*pii.^(k-1));
    end
end
gp = gs./ps*gsum;


%% Density
rho = 1./gp;
end