# BrineVolumeFraction
BrineVolumeFraction is a MATLAB code base for modeling the volume fraction of brine, salt, and ice (Ih) as a function of temperature and bulk ice salinity. Models are built by deriving phase behavior functions from the output of freezing simulations using the open source aqueous geochemistry software packages [FREZCHEM](https://github.com/MarcNeveu/frezchem) (v13.3, v15.1) and [PHREEQC](https://www.usgs.gov/software/phreeqc-version-3). 

## Related Publications
(3) Wolfenbarger, N. S., Blankenship, D. D., Young, D. A., Scanlan, K. M., Chivers, C. J., Findlay, D., Steinbruegge, G. B., Chan, K., Grima, C., Soderlund, K. M., & Schroeder, D. M. (XXXX). Radar Characterization of Salt Layers in Europa’s Ice Shell as a Window into Critical Ice-Ocean Exchange Processes.

(2) Wolfenbarger, N. S., Fox-Powell, M. G., Buffo, J. J., Soderlund, K. M., & Blankenship, D. D. (2022). Brine Volume Fraction as a Habitability Metric for Europa's Ice Shell. Geophysical Research Letters, 49(22), e2022GL100586, https://doi.org/10.1029/2022GL100586.

(1) Wolfenbarger, N. S., Fox-Powell, M. G., Buffo, J. J., Soderlund, K. M., & Blankenship, D. D. (2022). Compositional controls on the distribution of brine in Europa's ice shell. Journal of Geophysical Research: Planets, 127(9), e2022JE007305. https://doi.org/10.1029/2022JE007305 


## Functions
```
FREZCHEM\
[P,T,rho_b] = brine_density_FREZCHEM(fn)
[P,T,F1,F2,F3,k_star,rho_ss] = F1F2F3_FREZCHEM(fn)
[P,T,Sb] = liquidus_FREZCHEM(fn)  
FrOut = read_FrOut(fn)  
FrOut = read_FrOut133(fn)  
FrOut = read_FrOut151(fn) 
[P,Tmat,Smat,Vi_V, Vb_V, Vs_V] = volume_fraction_FREZCHEM(T,S,fn)  

Ice Ih Properties\
T = Tmelt(P)  
T = Tsub(P)
T = Tsat(P)
P = Pmelt(T)
P = Psub(T)
P = Psat(T)
rho = ice_density(T,P)  

PHREEQC\
[P,T,rho_b] = brine_density_PHREEQC(fn)
[P,T,F1,F2,F3,k_star,rho_ss] = F1F2F3_PHREEQC(fn)
[P,T,Sb] = liquidus_PHREEQC(fn)
pqo = read_pqo(fn)  
[P,Tmat,Smat,Vi_V, Vb_V, Vs_V] = volume_fraction_PHREEQC(T,S,fn)

Water Properties\
rho = water_density(T,P)  
``` 
## Function Descriptions
`brine_density_FREZCHEM(fn)` extracts the temperature and brine density (where ice is thermodynamically stable) from the output file of a FREZCHEM (v13.3,v15.1) or PHREEQC simulation, including a data point at the pure ice pressure melting temperature obtained from `Tmelt(P)`.

`brine_density_PHREEQC(fn)` extracts the pressure (Pa), temperature (C) and brine density (g/cm^3) for a solution composition and pressure defined by the PHREEQC output over the temperature range from the freezing point to the eutectic point.

`F1F2F3_FREZCHEM(fn)` derives phase behavior functions from the FREZCHEM (v13.3,v15.1) or PHREEQC simulation output file specified by `fn`.

`F1F2F3_PHREEQC(fn)` calculates the phase behavior functions `F1` and `F2` in Cox and Weeks (1983) and a new phase behavior function `F3` from freezing simulations in PHREEQC over a temperature range from the pure ice pressure melting point to the eutectic point.

`ice_density(T,P)` calculates the density of ice Ih using the equations of state, originally published by Wagner et al. (2011), included in the International Association for the Properties of Water and Steam (IAPWS) Revised Release on the Equation of State 2006 for H2O Ice Ih: IAPWS R10-06(2009).

`liquidus_FREZCHEM(fn)` extracts the temperature and brine salinity (where ice is stable) from the FREZCHEM (v13.3,v15.1) or PHREEQC simulation output file specified by `fn`, including a data point at the pure ice pressure melting temperature obtained from `Tmelt(P)`. 

`liquidus_PHREEQC(fn)` extracts the temperature (C) and brine salinity (ppt) which define the liquidus curve for a solution composition and pressure defined by PHREEQC over the temperature range from the freezing point to the eutectic point.

`Pmelt(T)` calculates the melting pressure `P` of ice Ih as a function of temperature `T` using the equations of state originally published by Wagner et al. (2011), included in the International Association for the Properties of Water and Steam (IAPWS) Revised Release on the Pressure along the Melting and Sublimation Curves of Ordinary Water Substance: IAPWS R14-08(2011).

`Psat(T)` calculates the saturation pressure `P` of water given temperature `T` using the equations included in the International Association for the Properties of Water and Steam (IAPWS) Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam: IAPWS R7-97(2012).

`Psub(T)` calculates the sublimation pressure `P` of ice Ih given temperature `T` using the equations of state originally published by Wagner et al. (2011), included in the International Association for the Properties of Water and Steam (IAPWS) Revised Release on the Pressure along the Melting and Sublimation Curves of Ordinary Water Substance: IAPWS R14-08(2011).

`read_FrOut(fn)`  reads a FREZCHEM output file by calling either `read_FrOut133(fn)` or `read_FrOut151(fn)`, depending on the version of the output file, to create the structure `FrOut`. 

`read_FrOut133(fn)` creates a structure `FrOut` from a FREZCHEM v13.3 output file `fn` that contains bulk properties, solution properties, and solid properties at each convergent temperature step. It also estimates solid properties at the eutectic temperature (as opposed to the final convergent temperature step).

`read_FrOut151(fn)` creates a structure `FrOut` from a FREZCHEM v15.1 output file `fn` that contains bulk properties, solution properties, and solid properties at each convergent temperature step. It also estimates solid properties at the eutectic temperature (as opposed to the final convergent temperature step).

`read_pqo(fn)` extracts data from PHREEQC output file, assuming input file decreases reaction temperature in steps.
 
`Tmelt(P)` calculates the melting temperature `T` of ice Ih as a function of pressure `P` using the equations originally published by Wagner et al. (2011), included in the International Association for the Properties of Water and Steam (IAPWS) Revised Release on the Pressure along the Melting and Sublimation Curves of Ordinary Water Substance: IAPWS R14-08(2011).
 
`Tsat(P)` calculates the saturation temperature `T` of water given pressure `P` using the equations included in the International Association for the Properties of Water and Steam (IAPWS) Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam: IAPWS R7-97(2012).
 
`Tsub(P)` calculates the sublimation temperature `T` of ice Ih given pressure `P` using the equations originally published by Wagner et al. (2011), included in the International Association for the Properties of Water and Steam (IAPWS) Revised Release on the Pressure along the Melting and Sublimation Curves of Ordinary Water Substance: IAPWS R14-08(2011).
 
`volume_fraction_FREZCHEM(T,S,fn)` calculates the volume fraction of ice `Vi_V`, brine `Vb_V`, and solid salt `Vs_V` as a function of temperature `T` and bulk ice salinity `S`, using the phase behavior functions derived from the FREZCHEM output file `fn`.

 
`volume_fraction_PHREEQC(T,S,fn)` calculates the volume fraction of ice `Vi_V`, brine `Vb_V`, and solid salt `Vs_V` as a function of temperature `T` and bulk ice salinity `S`, using the phase behavior functions derived from the PHREEQC output file `fn`.

`water_density(T,P)` calculates the density of water as a function of temperature `T` and pressure `P` using the equations of state, originally published by Feistel et al. (2003, 2008), included in the International Association for the Properties of Water and Steam (IAPWS) Supplementary Release on a Computationally Efficient Thermodynamic Formulation for Liquid Water for Oceanographic Use: IAPWS SR7-09. 

