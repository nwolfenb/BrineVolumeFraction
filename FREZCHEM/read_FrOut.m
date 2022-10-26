function FrOut = read_FrOut(fn)
% Extracts data from FREZCHEM output file (v13.3 or v15.1)
%
% Syntax: 
% FrOut = read_FrOut(fn)
%
% Inputs:
% fn        filename
%
% Outputs:
%
% FrOut     structure
%-----------Size-----------------------------------------------------------
%           N = number of temperature points
%           M = number of solution species
%           L = number of solid species
%
%-----------Bulk Properties (N x 1)----------------------------------------
%           FrOut.T                 Temperatures (K)
%           FrOut.ion_str           Ionic Strength
%           FrOut.rho_b             Brine Density (g/cm^3)
%           FrOut.osmotic_coef      Osmotic Coefficient
%           FrOut.br_H2O            Mass of Water in Brine (g)
%           FrOut.ice_H2O           Mass of Water in Ice (g)
%           FrOut.pressure          Pressure (bar)
%           FrOut.br_salt           Mass of Salt (Ions) in Brine (g)
%           FrOut.sm                Mass of Solid Salt (Hydrates) in Brine
%                                   (g)
%           FrOut.rho_ss            Mean Density of Solid Salt (Hydrates)
%                                   (g/cm^3)
%
%-----------Solution Properties (N x M)------------------------------------
%           FrOut.solution_species  Solution Species (Ions)
%           FrOut.init_conc         Initial Concentration (moles/kg)
%           FrOut.final_conc        Final Concentration(moles/kg)
%           FrOut.activity_coef     Activity Coefficient
%           FrOut.activity          Activity
%           FrOut.solution_moles	Moles of Solution Species (moles)
%           FrOut.mass_balance      Mass of Major Species (moles/kg)
%
%-----------Solid Properties (N x L)---------------------------------------
%           FrOut.solid_species     Solid Species (Hydrates)
%           FrOut.solid_moles       Solid Moles (moles)
%           FrOut.equil_const       Equilibrium Constant
%           FrOut.accum_moles       Accumulated Moles (moles)
%
%-----------Salts Precipitating at the Eutectic (L* x 1)-------------------
% NOTE: This field will be populated with NaN if FREZCHEM does not reach
% the eutectic (sometimes it doesn't due to convergence issues). 
%
%           FrOut.eutectic.T                Eutectic Temperature*
%           FrOut.eutectic.solid_species    Solid Species (Hydrates) 
%                                           Precipitating at the Eutectic
%           FrOut.eutectic.fw               Weight Fraction Solid Salt
%                                           (Hydrates) Precipitating at 
%                                           the Eutectic
%           FrOut.eutectic.rho_ss           Mean Density of Solid Salt
%                                           (Hydrates) Precipitating at 
%                                           the Eutectic
%
% *If the eutectic temperature is defined as the lowest temperature at
% which melt can occur, this temperature corresponds to the temperature one
% step below that (i.e, where no melt is stable).
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%% Detect FREZCHEM version from FREZCHEM output file
% Find number of solid species
flag = false;
k_solid = 0;

fid = fopen(fn);
tline = fgetl(fid);
while ischar(tline)
    tline = fgetl(fid);
    if ~isempty(tline) & tline~=-1
        if contains(tline,' Solid')
            tline = fgetl(fid);
            while ~isempty(tline)
                tline = fgetl(fid);
                k_solid = k_solid + 1;
            end
            flag = true;
        end
    end
    if flag
        break
    end
end
fclose(fid);
k_solid_tot = k_solid - 1;

disp(['File: ',fn])
if k_solid_tot == 101
    disp('FREZCHEM v13.3')
    FrOut = read_FrOut133(fn);
elseif k_solid_tot == 108
    disp('FREZCHEM v15.1')
    FrOut = read_FrOut151(fn);
else
    error('FREZCHEM version not recognized. readFrOut is only compatible with FREZCHEM versions 13.3 and 15.1.')
end

end

