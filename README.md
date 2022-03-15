# BrineVolumeFraction
MATLAB code base for modeling the volume fraction of brine, salt, and ice (Ih) as a function of temperature and bulk ice salinity using phase behavior functions derived from the output of freezing simulations using FREZCHEM v13.3 and 15.1.

The code base includes 7 functions:
FrOut = read_FrOut(fn)
FrOut = read_FrOut133(fn)
FrOut = read_FrOut151(fn)
[T,Sb] = liquidus_FREZCHEM(fn)
[Tmat,Smat,Vi_V, Vb_V, Vs_V] = volume_fraction_FREZCHEM(T,S,fn)

FrOut = read_FrOut(fn)
read_FrOut reads a FREZCHEM output file and calls either read_FrOut133(fn) or read_FrOut151(fn), depending on the version of the output file, to create the structure FrOut. 

FrOut = read_FrOut133(fn)
read_FrOut133(fn) creates a structure (FrOut) from a FREZCHEM v13.3 output file (fn) that contains bulk properties, solution properties, and solid properties at each convergent temperature step. It also estimates solid properties at the eutectic temperature (as opposed to the final convergent temperature step).

FrOut = read_FrOut151(fn)
read_FrOut151(fn) creates a structure (FrOut) from a FREZCHEM v15.1 output file (fn) that contains bulk properties, solution properties, and solid properties at each convergent temperature step. It also estimates solid properties at the eutectic temperature (as opposed to the final convergent temperature step).

[T,Sb] = liquidus_FREZCHEM(fn)
liquidus_FREZCHEM extracts the temperature and brine salinity (where ice is stable) from the output file of a FREZCHEM v13.3 or 15.1 simulation, including a data point at the pure ice pressure melting temperature obtained from Tmelt(P). 

T = Tmelt(P)
Tmelt(P) estimates the pressure melting temperature of ice Ih using the equations, originally published by Wagner et al. (2011), included in the International Association for the Properties of Water and Steam (IAPWS) Revised Release on the Equation of State 2006 for H2O Ice Ih: IAPWS R10-06(2009).

