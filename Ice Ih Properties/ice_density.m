function rho = ice_density(T,P)
% Density of ice Ih as a function of temperature and pressure from the
% Gibbs energy equation of state. Densities where ice is not
% thermodynamically stable will putput as NaN.
%
% Syntax:
% rho = ice_density(T,P)
%
% Inputs:
% T     Temperature (K), scalar or vector
% P     Pressure (Pa), scalar or vector
%
% Outputs:
% rho   Density (kg/m^3), scalar or vector
%
% Source:
% IAPWS R10-06(2009)
% http://www.iapws.org/relguide/Ice-Rev2009.pdf
%
% Range of Validity:
% 0 K <= T <= 273.16 K
% 0 < P <= 210 MPa
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%% Column Vectors
flip = false;
if ~iscolumn(T) && ~ismatrix(T)
    T = T';
    flip = true;
end

if ~iscolumn(P) && ~ismatrix(P)
    P = P';
    flip = true;
end

if length(P)~= length(T) && (~isscalar(P) && ~isscalar(T))
    error('Pressure/temperature input must either be a vector the same length as the temperature/pressure or a scalar value that applies to all temperatures/pressures.')
end

%% Coefficients of the equation of state (Gibbs potential function)
g00 = -0.632020233335886e6; % J/kg
g01 = 0.655022213658955; % J/kg
g02 = -0.189369929326131e-7; % J/kg
g03 = 0.339746123271053e-14; % J/kg
g04 = -0.556464869058991e-21; % J/kg
g0k = [g00 g01 g02 g03 g04];

t2 = 0.337315741065416+1i*0.335449415919309; % J/kg/K

r20 = -0.725974574329220e2+1i*-0.781008427112870e2; % J/kg/K
r21 = -0.557107698030123e-4+1i*0.464578634580806e-4; % J/kg/K
r22 = 0.234801409215913e-10+1i*-0.285651142904972e-10; % J/kg/K
r2k = [r20 r21 r22];

%% Constants
Tt = 273.16; % K
Pt = 611.657; % Pa
pii0 = 101325/Pt;

%% Dimensionless variables
tau = T/Tt;
pii = P/Pt;

%% g0p
g0p = 0;
for k = 1:4
    g0p = g0p + g0k(k+1)*(k/Pt)*(pii-pii0).^(k-1);
end

%% r2p
r2p = 0;
for k = 1:2
    r2p = r2p + r2k(k+1)*(k/Pt)*(pii-pii0).^(k-1);
end

%% Derivative of Gibbs energy, gp(T,P)
Re = real(r2p.*((t2-tau).*log(t2-tau)+(t2+tau).*log(t2+tau)-2*t2*log(t2)-tau.^2/t2));
gp = g0p + Tt.*Re;

%% Density
rho = 1./gp;

if flip
    rho = rho';
end

%% Check phase stability 
T_ref = NaN(size(P));
T_ref(P>Pt) = Tmelt(P(P>Pt));
T_ref(P<Pt) = Tsub(P(P<Pt));
T_ref(P==Pt) = Tt;
rho(T>T_ref) = NaN;
rho(isnan(T_ref)) = NaN;

if any(isnan(rho))
    warning('Ice is not thermodyanmically stable for at least one input temperature and pressure.')
end

end