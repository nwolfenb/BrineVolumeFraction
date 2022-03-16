# BrineVolumeFraction
MATLAB code base for modeling the volume fraction of brine, salt, and ice (Ih) as a function of temperature and bulk ice salinity using phase behavior functions derived from the output of freezing simulations using FREZCHEM v13.3 and 15.1. This code base is used to generate the results in the manuscript "Compositional Controls on the Distribution of Brine in Europa's Ice Shell", submitted to JGR: Planets.

The code base includes the following functions:  
T = Tmelt(p)  
rho = water_density(T,p)  
rho = ice_density(T,p)  
FrOut = read_FrOut(fn)  
FrOut = read_FrOut133(fn)  
FrOut = read_FrOut151(fn)  
[T,Sb] = liquidus_FREZCHEM(fn)  
[Tmat,Smat,Vi_V, Vb_V, Vs_V] = volume_fraction_FREZCHEM(T,S,fn)  

T = Tmelt(p)  
Tmelt(p) calculates the melting temperature of ice as a function of pressure (p) Ih using the equations of state, originally published by Wagner et al. (2011), included in the International Association for the Properties of Water and Steam (IAPWS) Revised Release on the Equation of State 2006 for H2O Ice Ih: IAPWS R10-06(2009).

rho = water_density(T,p)  
water_density(T,p) calculates the density of water as a function of temperature (T) and pressure (p) using the equations of state, originally published by Feistel et al. (2003, 2008), included in the International Association for the Properties of Water and Steam (IAPWS) Supplementary Release on a Computationally Efficient Thermodynamic Formulation for Liquid Water for Oceanographic Use: IAPWS SR7-09. 

rho = ice_density(T,p)  
ice_density(T,p) calculates the density of ice Ih using the equations of state, originally published by Wagner et al. (2011), included in the International Association for the Properties of Water and Steam (IAPWS) Revised Release on the Equation of State 2006 for H2O Ice Ih: IAPWS R10-06(2009).

FrOut = read_FrOut(fn)  
read_FrOut reads a FREZCHEM output file and calls either read_FrOut133(fn) or read_FrOut151(fn), depending on the version of the output file, to create the structure FrOut. 

FrOut = read_FrOut133(fn)  
read_FrOut133(fn) creates a structure (FrOut) from a FREZCHEM v13.3 output file (fn) that contains bulk properties, solution properties, and solid properties at each convergent temperature step. It also estimates solid properties at the eutectic temperature (as opposed to the final convergent temperature step).

FrOut = read_FrOut151(fn)  
read_FrOut151(fn) creates a structure (FrOut) from a FREZCHEM v15.1 output file (fn) that contains bulk properties, solution properties, and solid properties at each convergent temperature step. It also estimates solid properties at the eutectic temperature (as opposed to the final convergent temperature step).

[T,Sb] = liquidus_FREZCHEM(fn)  
liquidus_FREZCHEM extracts the temperature and brine salinity (where ice is stable) from the output file of a FREZCHEM v13.3 or 15.1 simulation, including a data point at the pure ice pressure melting temperature obtained from Tmelt(p). 

[Tmat,Smat,Vi_V, Vb_V, Vs_V] = volume_fraction_FREZCHEM(T,S,fn)  
volume_fraction_FREZCHEM(T,S,fn) calculates the volume fraction of ice (Vi_V), brine (Vb_V), and solid salt (Vs_V) as a function of temperature (T) and bulk ice salinity (S), using the phase behavior functions derived from the FREZCHEM output file (fn).
