function pqo = read_pqo(fn)
% Extracts data from PHREEQC output file. File format assumes reaction
% temperature is decreased using
%
% REACTION_TEMPERATURE
%        [Ti] [Tf] in [X] steps
%
% in PHREEQC input file.
%
% Syntax:
% pqo = read_pqo(fn)
%
% Inputs:
% fn        Filename of PHREEQC output file, string
%
% Outputs:
% pqo    structure
%-----------Size-----------------------------------------------------------
%           N = number of temperature steps
%           M = number of elements
%           L = number of solid phases
%           P = number of solution species
%
%-----------Derived bulk properties (N x 1)--------------------------------
%           pqo.br_H2O              Mass of Water in Brine (g)
%           pqo.ice_H2O             Mass of Water in Ice (g)
%           pqo.br_salt             Mass of Salt (Ions) in Brine (g)
%           pqo.sm                  Mass of Solid Salt (Hydrates) in Brine
%                                   (g)
%           pqo.rho_ss              Mean Density of Solid Salt (Hydrates)
%                                   (g/cm^3)
%
%-----------Distribution of species (N x P)--------------------------------
%           pqo.species.name        Name of species
%           pqo.species.molal       Molality of species (moles/kgw)
%
%-----------Description of solution (N x 1)--------------------------------
%           pqo.pH                  pH
%           pqo.pe                  pe
%           pqo.cond                Specific Conductance (uS/cm)
%           pqo.rho_b               Brine Density (g/cm^3)
%           pqo.vol_H2O             Volume of water (L)
%           pqo.activity_H2O        Activity of water
%           pqo.ion_str             Ionic Strength (mol/kgw)
%           pqo.mass_H2O            Mass of water (kg)
%           pqo.alk                 Total alkalinity (eq/kg)
%           pqo.T                   Temperature (K)
%           pqo.pressure            Pressure (atm)
%           pqo.E_balance           Electrical balance (eq)
%           pqo.error               Percent error, 100*(Cat-|An|)/(Cat+|An|)
%           pqo.iterations          Iterations
%           pqo.gamma_iterations    Gamma Iterations
%           pqo.osmotic_coef        Osmotic coefficient
%           pqo.rho_H2O             Density of water (g/cm^3)
%
%-----------Solution composition (N x M)-----------------------------------
%           pqo.elements.name       Name of Elements
%           pqo.elements.molal      Molality of Elements (moles/kgw)
%           pqo.elements.moles      Moles of Elements (moles)
%
%-----------Phase assemblage (N x L)---------------------------------------
%           pqo.solids            Name of Solid Phase
%           pqo.SI                Saturation Index
%           pqo.logIAP            Ion Activity Product
%           pqo.logK              Equilibrium Constant
%           pqo.solid_moles       Moles in Assemblage (Final)
%
%-----------Salts precipitating at the eutectic (L* x 1)-------------------
%           pqo.eutectic.T                  Eutectic Temperature*
%           pqo.eutectic.solid_species      Solid Species (Hydrates)
%                                           Precipitating at the Eutectic
%           pqo.eutectic.fw                 Weight Fraction Solid Salt
%                                           (Hydrates) Precipitating at
%                                           the Eutectic
%           pqo.eutectic.rho_ss             Mean Density of Solid Salt
%                                           (Hydrates) Precipitating at
%                                           the Eutectic
%
% *For now, the eutectic temperature is defined as the temperature where
% the mass of brine is less than 1e-3 kg
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%% Elemental Properties (PHREEQC)

% SOLUTION_MASTER_SPECIES
symb = {'H';...
    'O';...
    'Ca';...
    'Mg';...
    'Na';...
    'K';...
    'Cl';...
    'C';...
    'S'};

mass = [1.008;... % H
    15.999;... % O
    40.08;... % Ca
    24.31;... % Mg
    22.99;... % Na
    39.1;... % K
    35.45;... % Cl
    12.015;... % C
    32.064]; % S

% Species Name
species.str = {'H+';...
    'H2O';...
    'Ca+2';...
    'Mg+2';...
    'Na+';...
    'K+';...
    'Cl-';...
    'CO3-2';...
    'SO4-2';...
    'OH-';...
    'HCO3-';...
    'CO2';...
    'MgOH+'};

% Species Mass (g)
species.mass = [mass(strcmp('H',symb));... % H+
    2*mass(strcmp('H',symb))+mass(strcmp('O',symb));... % H2O
    mass(strcmp('Ca',symb));... % Ca+2
    mass(strcmp('Mg',symb)); ... % Mg+2
    mass(strcmp('Na',symb)); ... % Na+
    mass(strcmp('K',symb)); ... % K+
    mass(strcmp('Cl',symb)); ... % Cl-
    mass(strcmp('C',symb))+3*mass(strcmp('O',symb)); ... % CO3-2
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb)); ... % SO4-2
    mass(strcmp('O',symb))+mass(strcmp('H',symb)); ... % OH-
    mass(strcmp('H',symb))+mass(strcmp('C',symb))+3*mass(strcmp('O',symb)); ... % HCO3-
    mass(strcmp('C',symb))+2*mass(strcmp('O',symb)); ... % CO2
    mass(strcmp('Mg',symb))+mass(strcmp('O',symb))+mass(strcmp('H',symb))]; ... % MgOH
    
solids.str = {'Anhydrite';... % CaSO4
    'Antarcticite';... % CaCl2:6H2O
    'Aphthitalite';...  % Na2SO4:3K2SO4
    'Aragonite';... % CaCO3
    'Arcanite';... % K2SO4
    'Bassanite';... % CaSO4:0.5H2O
    'Bischofite';... % MgCl2:6H2O
    'Bloedite';... % Na2Mg(SO4)2:4H2O
    'Calcite';... % CaCO3
    'Carnallite';... % KMgCl3:6H2O
    'Dolomite';... % CaMg(CO3)2
    'Epsomite';... % MgSO4:7H2O
    'Glaserite';... % Na2SO4:3K2SO4
    'Glauberite';... % Na2SO4:CaSO4
    'Gypsum';... % CaSO4:2H2O
    'Halite';... % NaCl
    'Hexahydrite';... % MgSO4:6H2O
    'Hydrohalite';... % NaCl:2H2O
    'Hydromagnesite';... % (MgCO3)3:Mg(OH)2:3H2O
    'Ice(s)';... % H2O
    'Ikaite';... % CaCO3:6H2O
    'Kainite';... % KCl:MgSO4:2.75H2O
    'Kalicinite';... % KHCO3
    'Kieserite';... % MgSO4:H2O
    'Labile_Salt';... % Na4(SO4)2:CaSO4:2H2O
    'Landsfordite';... % MgCO3:5H2O
    'Langbeinite';... % K2SO4:2MgSO4
    'Leonite';... % K2SO4:MgSO4:4H2O
    'Magnesite';... % MgCO3
    'Meridianite';... % MgSO4:11H2O
    'Meridianiite';... % MgSO4:11H2O
    'MgCl2:12H2O';... % MgCl2:12H2O
    'MgCl2:8H2O';... % MgCl2:8H2O
    'Mirabilite';... % Na2SO4:10H2O
    'Na2CO3:7H2O';... % Na2CO3:7H2O
    'Na2SO4:7H2O';... % Na2SO4:7H2O
    'Nahcolite';... % NaHCO3
    'Natron';... % Na2CO3:10H2O
    'Nesquehonite';... % MgCO3:3H2O
    'Pentahydrite';... % MgSO4:5H2O
    'Picromerite';... % MgSO4:K2SO4:6H2O
    'Polyhalite';... % K2SO4:MgSO4:Ca2(SO4)2:2H2O
    'Starkeyite';... % MgSO4:4H2O
    'Sylvite';... % KCl
    'Syngenite';... % K2SO4:CaSO4:H2O
    'Tachyhydrite';... % CaCl2:(MgCl2)2:12H2O
    'Thenardite';... % Na2SO4
    'Trona';... % Na3H(CO3)2:2H2O
    'Vaterite'}; % CaCO3

solids.mass = [mass(strcmp('Ca',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb));... % CaSO4
    mass(strcmp('Ca',symb))+2*mass(strcmp('Cl',symb))+...
    6*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % CaCl2:6H2O
    2*mass(strcmp('Na',symb))+6*mass(strcmp('K',symb))+...
    4*(mass(strcmp('S',symb))+4*mass(strcmp('O',symb)));...  % Na2SO4:3K2SO4
    mass(strcmp('Ca',symb))+...
    mass(strcmp('C',symb))+3*mass(strcmp('O',symb));... % CaCO3
    2*mass(strcmp('K',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb));... % K2SO4
    mass(strcmp('Ca',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    0.5*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % CaSO4:0.5H2O
    mass(strcmp('Mg',symb))+2*mass(strcmp('Cl',symb))+...
    6*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % MgCl2:6H2O
    2*mass(strcmp('Na',symb))+mass(strcmp('Mg',symb))+...
    2*(mass(strcmp('S',symb))+4*mass(strcmp('O',symb)))+...
    4*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % Na2Mg(SO4)2:4H2O
    mass(strcmp('Ca',symb))+...
    mass(strcmp('C',symb))+3*mass(strcmp('O',symb));... % CaCO3
    mass(strcmp('K',symb))+mass(strcmp('Mg',symb))+...
    3*mass(strcmp('Cl',symb))+...
    6*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % KMgCl3:6H2O
    mass(strcmp('Ca',symb))+mass(strcmp('Mg',symb))+...
    2*(mass(strcmp('C',symb))+3*mass(strcmp('O',symb)));... % CaMg(CO3)2
    mass(strcmp('Mg',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    7*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % MgSO4:7H2O
    2*mass(strcmp('Na',symb))+6*mass(strcmp('K',symb))+...
    4*(mass(strcmp('S',symb))+4*mass(strcmp('O',symb)));... % Na2SO4:3K2SO4
    2*mass(strcmp('Na',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    mass(strcmp('Ca',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb));... % Na2SO4:CaSO4
    mass(strcmp('Ca',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    2*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % CaSO4:2H2O
    mass(strcmp('Na',symb))+mass(strcmp('Cl',symb));... % NaCl
    mass(strcmp('Mg',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    6*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % MgSO4:6H2O
    mass(strcmp('Na',symb))+mass(strcmp('Cl',symb))+...
    2*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % NaCl:2H2O
    3*(mass(strcmp('Mg',symb))+...
    mass(strcmp('C',symb))+3*mass(strcmp('O',symb)))+...
    mass(strcmp('Mg',symb))+...
    2*(mass(strcmp('O',symb))+mass(strcmp('H',symb)))+...
    3*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % (MgCO3)3:Mg(OH)2:3H2O
    2*mass(strcmp('H',symb))+mass(strcmp('O',symb));... % H2O
    mass(strcmp('Ca',symb))+...
    mass(strcmp('C',symb))+3*mass(strcmp('O',symb))+...
    6*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % CaCO3:6H2O
    mass(strcmp('K',symb))+mass(strcmp('Cl',symb))+...
    mass(strcmp('Mg',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    2.75*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % KCl:MgSO4:2.75H2O
    mass(strcmp('K',symb))+mass(strcmp('H',symb))+...
    mass(strcmp('C',symb))+3*mass(strcmp('O',symb));... % KHCO3
    mass(strcmp('Mg',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    2*mass(strcmp('H',symb))+mass(strcmp('O',symb));... % MgSO4:H2O
    4*mass(strcmp('Na',symb))+...
    2*(mass(strcmp('S',symb))+4*mass(strcmp('O',symb)))+...
    mass(strcmp('Ca',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    2*mass(strcmp('H',symb))+mass(strcmp('O',symb));... % Na4(SO4)2:CaSO4:2H2O
    mass(strcmp('Mg',symb))+...
    mass(strcmp('C',symb))+3*mass(strcmp('O',symb))+...
    5*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % MgCO3:5H2O
    2*mass(strcmp('K',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    2*(mass(strcmp('Mg',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb)));... % K2SO4:2MgSO4
    2*mass(strcmp('K',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    mass(strcmp('Mg',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    4*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % K2SO4:MgSO4:4H2O
    mass(strcmp('Mg',symb))+...
    mass(strcmp('C',symb))+3*mass(strcmp('O',symb));... % MgCO3
    mass(strcmp('Mg',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    11*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % MgSO4:11H2O
    mass(strcmp('Mg',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    11*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % MgSO4:11H2O
    mass(strcmp('Mg',symb))+2*mass(strcmp('Cl',symb))+...
    12*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % MgCl2:12H2O
    mass(strcmp('Mg',symb))+2*mass(strcmp('Cl',symb))+...
    8*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % MgCl2:8H2O
    2*mass(strcmp('Na',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    10*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % Na2SO4:10H2O
    2*mass(strcmp('Na',symb))+...
    mass(strcmp('C',symb))+3*mass(strcmp('O',symb))+...
    7*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % Na2CO3:7H2O
    2*mass(strcmp('Na',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    7*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % Na2SO4:7H2O
    mass(strcmp('Na',symb))+mass(strcmp('H',symb))+...
    mass(strcmp('C',symb))+3*mass(strcmp('O',symb));... % NaHCO3
    2*mass(strcmp('Na',symb))+...
    mass(strcmp('C',symb))+3*mass(strcmp('O',symb))+...
    10*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % Na2CO3:10H2O
    mass(strcmp('Mg',symb))+...
    mass(strcmp('C',symb))+3*mass(strcmp('O',symb))+...
    3*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % MgCO3:3H2O
    mass(strcmp('Mg',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    5*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % MgSO4:5H2O
    mass(strcmp('Mg',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    2*mass(strcmp('K',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    6*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % MgSO4:K2SO4:6H2O
    2*mass(strcmp('K',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    mass(strcmp('Mg',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    2*mass(strcmp('Ca',symb))+...
    2*(mass(strcmp('S',symb))+4*mass(strcmp('O',symb)))+...
    2*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % K2SO4:MgSO4:Ca2(SO4)2:2H2O
    mass(strcmp('Mg',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    4*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % MgSO4:4H2O
    mass(strcmp('K',symb))+mass(strcmp('Cl',symb));... % KCl
    2*mass(strcmp('K',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    mass(strcmp('Ca',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb))+...
    2*mass(strcmp('H',symb))+mass(strcmp('O',symb));... % K2SO4:CaSO4:H2O
    mass(strcmp('Ca',symb))+2*mass(strcmp('Cl',symb))+...
    2*(mass(strcmp('Mg',symb))+2*mass(strcmp('Cl',symb)))+...
    12*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % CaCl2:(MgCl2)2:12H2O
    2*mass(strcmp('Na',symb))+...
    mass(strcmp('S',symb))+4*mass(strcmp('O',symb));... % Na2SO4
    3*mass(strcmp('Na',symb))+mass(strcmp('H',symb))+...
    2*(mass(strcmp('C',symb))+3*mass(strcmp('O',symb)))+...
    2*(2*mass(strcmp('H',symb))+mass(strcmp('O',symb)));... % Na3H(CO3)2:2H2O
    mass(strcmp('Ca',symb))+...
    mass(strcmp('C',symb))+3*mass(strcmp('O',symb))]; % CaCO3

% Solids Molar Volume (cm^3/mol)
solids.molvol =  [45.94;... % CaSO4
    128.12;... % CaCl2:6H2O
    246.24;...  % Na2SO4:3K2SO4
    34.15;... % CaCO3
    65.50;... % K2SO4
    53.03;... % CaSO4:0.5H2O
    129.57;... % MgCl2:6H2O
    149.98;... % Na2Mg(SO4)2:4H2O
    36.93;... % CaCO3
    172.58;... % KMgCl3:6H2O
    64.34;... % CaMg(CO3)2
    146.71;... % MgSO4:7H2O
    246.24;... % Na2SO4:3K2SO4
    100.82;... % Na2SO4:CaSO4
    74.69;... % CaSO4:2H2O
    27.02;... % NaCl
    132.58;... % MgSO4:6H2O
    57.96;... % NaCl:2H2O
    169.13;... % (MgCO3)3:Mg(OH)2:3H2O
    19.65;... % H2O
    117.54;... % CaCO3:6H2O
    114.40;... % KCl:MgSO4:2.75H2O
    46.14;... % KHCO3
    56.60;... % MgSO4:H2O
    NaN;... % Na4(SO4)2:CaSO4:2H2O
    100.80;... % MgCO3:5H2O
    146.6431;... % K2SO4:2MgSO4
    166.6015;... % K2SO4:MgSO4:4H2O
    28.02;... % MgCO3
    207.44;... % MgSO4:11H2O
    207.44;... % MgSO4:11H2O
    218.10;... % MgCl2:12H2O
    159.08;... % MgCl2:8H2O
    219.80;... % Na2SO4:10H2O
    153.71;... % Na2CO3:7H2O
    171.7803;... % Na2SO4:7H2O
    38.91;... % NaHCO3
    198.71;... % Na2CO3:10H2O
    74.79;... % MgCO3:3H2O
    110.9942;... % MgSO4:5H2O
    191.78;... % MgSO4:K2SO4:6H2O
    218.27;... % K2SO4:MgSO4:Ca2(SO4)2:2H2O
    95.87;... % MgSO4:4H2O
    37.52;... % KCl
    127.8;... % K2SO4:CaSO4:H2O
    311.81;... % CaCl2:(MgCl2)2:12H2O
    53.33;... % Na2SO4
    107.02;... % Na3H(CO3)2:2H2O
    37.72]; % CaCO3

% Species Density (g/cm^3)
solids.rho = solids.mass./solids.molvol;

%% Read-In File
% Find number of temperature steps
fid = fopen(fn);
tline = fgetl(fid);
N = 0;
flag = 0;
while ischar(tline)
    tline = fgetl(fid);
    if ~isempty(tline) & tline~=-1
        if contains(tline,'Reaction step')
            while ~contains(tline,'-Phase assemblage-')
                if contains(tline,'ERROR')
                    flag = 1;
                end
                tline = fgetl(fid);
            end
            if flag == 1
                break
            end
            N = N+1;
            fgetl(fid);
            
            tline = fgetl(fid);
            L = 0;
            while ~isempty(tline)
                L = L+1;
                tline = fgetl(fid);
            end
            
            while ~contains(tline,'-Solution composition-')
                tline = fgetl(fid);
            end
            for i = 1:3
                fgetl(fid);
            end
            tline = fgetl(fid);
            M = 0;
            while ~isempty(tline)
                M = M+1;
                tline = fgetl(fid);
            end
            
            while ~contains(tline,'-Distribution of species-')
                tline = fgetl(fid);
            end
            for i = 1:5
                fgetl(fid);
            end
            tline = fgetl(fid);
            while isspace(tline(1)) % Skip OH-, H+, and H2O
                tline = fgetl(fid);
            end
            P = 0;
            while ~isempty(tline)
                if isspace(tline(1))
                    P = P+1;
                end
                tline = fgetl(fid);
            end
        end
    end
end
fclose(fid);

% Pre-allocate

% Derived bulk properties
pqo.br_H2O =  NaN(N,1); % Mass of Water in Brine (g)
pqo.ice_H2O = pqo.br_H2O; % Mass of Water in Ice (g)
pqo.br_salt = pqo.br_H2O; % Mass of Salt (Ions) in Brine (g)
pqo.sm = pqo.br_H2O; % Mass of Solid Salt (Hydrates) in Brine (g)
pqo.rho_ss = pqo.br_H2O; % Mean Density of Solid Salt (Hydrates) (g/cm^3)

% Distribution of species
pqo.species.name = cell(N,P);
pqo.species.name(:) = {''}; % Name of species
pqo.species.molal = NaN(N,P); % Molality of species (moles/kgw)

% Description of solution
pqo.pH = NaN(N,1); % pH
pqo.pe = pqo.pH; % pe
pqo.cond = pqo.pH; % Specific Conductance (uS/cm)
pqo.rho_b = pqo.pH; % Brine Density (g/cm^3)
pqo.vol_H2O = pqo.pH; % Volume of water (L)
pqo.activity_H2O = pqo.pH; % Activity of water
pqo.ion_str = pqo.pH; % Ionic Strength (mol/kgw)
pqo.mass_H2O = pqo.pH; % Mass of water (kg)
pqo.alk = pqo.pH; % Total alkalinity (eq/kg)
pqo.T = pqo.pH; % Temperature (K)
pqo.pressure = pqo.pH; % Pressure (atm)
pqo.E_balance = pqo.pH; % Electrical balance (eq)
pqo.error = pqo.pH; % Percent error, 100*(Cat-|An|)/(Cat+|An|)
pqo.iterations = pqo.pH; % Iterations
pqo.gamma_iterations = pqo.pH; % Gamma Iterations
pqo.osmotic_coef = pqo.pH; % Osmotic coefficient
pqo.rho_H2O = pqo.pH; % Density of water (g/cm^3)

% Solution composition
pqo.elements.name = cell(N,M); % Name of Elements
pqo.elements.name(:) = {''};
pqo.elements.molal = zeros(N,M); % Molality (moles/kgw)
pqo.elements.moles = pqo.elements.molal; % Moles (moles)

% Phase assemblage
pqo.solids = cell(N,L); % Name of Solid Phase
pqo.solids(:) = {''};
pqo.SI = zeros(N,L); % Saturation Index
pqo.logIAP = pqo.SI; % Ion Activity Product
pqo.logK = pqo.SI; % Equilibrium Constant
pqo.solid_moles = pqo.SI; % Moles in Assemblage (Final)

% Collect Data
fid = fopen(fn);
tline = fgetl(fid);
n = 0;
flag = 0;
while ischar(tline)
    tline = fgetl(fid);
    if ~isempty(tline) & tline~=-1
        if contains(tline,'Reaction step')
            n = n+1;
            l = 0;
            m = 0;
            tline = fgetl(fid);
            while ~contains(tline,'-Phase assemblage-')
                if contains(tline,'ERROR')
                    flag = 1;
                end
                tline = fgetl(fid);
            end
            for i = 1:4
                fgetl(fid);
            end
            
            if flag == 1
                break
            end
            
            tline = fgetl(fid);
            while ~isempty(tline)
                l = l+1;
                solid_str = tline(1:14); % HERE
                pqo.solids{n,l} = solid_str(~isspace(solid_str));
                data = sscanf(tline,'%*s %f %f %f %*f %f %*f');
                if isempty(data)
                    pqo.SI(n,l) = NaN;
                    pqo.logIAP(n,l) = NaN;
                    pqo.logK(n,l) = NaN;
                    pqo.solid_moles(n,l) = NaN;
                else
                    pqo.SI(n,l) = data(1);
                    pqo.logIAP(n,l) = data(2);
                    pqo.logK(n,l) = data(3);
                    pqo.solid_moles(n,l) = data(4);
                end
                tline = fgetl(fid);
            end
            
            while ~contains(tline,'-Solution composition-')
                tline = fgetl(fid);
            end
            for i = 1:3
                fgetl(fid);
            end
            tline = fgetl(fid);
            while ~isempty(tline)
                m = m+1;
                solut_str = tline(1:3);
                pqo.elements.name{n,m} = solut_str(~isspace(solut_str));
                data = sscanf(tline,'%*s %f %f');
                pqo.elements.molal(n,m) = data(1);
                pqo.elements.moles(n,m) = data(2);
                tline = fgetl(fid);
            end
            
            while ~contains(tline,'-Description of solution-')
                tline = fgetl(fid);
            end
            fgetl(fid);
            tline = fgetl(fid);
            while ~isempty(tline)
                if contains(tline,'pH')
                    data = sscanf(tline,'%*s = %f %*s %*s');
                    pqo.pH(n) = data;
                elseif contains(tline,' pe ')
                    data = sscanf(tline,'%*s = %f %*s %*s %*s %*s');
                    pqo.pe(n) = data;
                elseif contains(tline,'Conductance')
                    data = sscanf(tline,'%*s %*s %*s %*s = %f');
                    pqo.cond(n) = data;
                elseif contains(tline,'Density') && ~contains(tline,'water')
                    data = sscanf(tline,'%*s %*s = %f');
                    pqo.rho_b(n) = data;
                elseif contains(tline,'Volume')
                    data = sscanf(tline,'%*s %*s = %f');
                    pqo.vol_H2O(n) = data;
                elseif contains(tline,'Activity')
                    data = sscanf(tline,'%*s %*s %*s = %f');
                    pqo.activity_H2O(n) = data;
                elseif contains(tline,'Ionic')
                    data = sscanf(tline,'%*s %*s %*s = %f');
                    pqo.ion_str(n) = data;
                elseif contains(tline,'Mass')
                    data = sscanf(tline,'%*s %*s %*s %*s = %f');
                    pqo.mass_H2O(n) = data;
                elseif contains(tline,'alkalinity')
                    data = sscanf(tline,'%*s %*s %*s = %f');
                    pqo.alk(n) = data;
                elseif contains(tline,'Temperature')
                    data = sscanf(tline,'%*s %*s = %f');
                    pqo.T(n) = data;
                elseif contains(tline,'Pressure')
                    data = sscanf(tline,'%*s %*s = %f');
                    pqo.pressure(n) = data;
                elseif contains(tline,'Electrical')
                    data = sscanf(tline,'%*s %*s %*s = %f');
                    pqo.E_balance(n) = data;
                elseif contains(tline,'Percent error')
                    data = sscanf(tline,'%*s %*s %*s = %f');
                    pqo.error(n) = data;
                elseif contains(tline,'Iterations') && ~contains(tline,'Gamma')
                    data = sscanf(tline,'%*s = %f %*s %*s');
                    pqo.iterations(n) = data;
                elseif contains(tline,'Gamma iterations')
                    data = sscanf(tline,'%*s %*s = %f');
                    pqo.gamma_iterations(n) = data;
                elseif contains(tline,'Osmotic')
                    data = sscanf(tline,'%*s %*s = %f');
                    pqo.osmotic_coef(n) = data;
                elseif contains(tline,'Density') && contains(tline,'water')
                    data = sscanf(tline,'%*s %*s %*s = %f');
                    pqo.rho_H2O(n) = data;
                end
                tline = fgetl(fid);
            end
            while ~contains(tline,'-Distribution of species-')
                tline = fgetl(fid);
            end
            for i = 1:5
                fgetl(fid);
            end
            tline = fgetl(fid);
            while isspace(tline(1)) % Skip OH-, H+, and H2O
                tline = fgetl(fid);
            end
            p = 0;
            while ~isempty(tline)
                if isspace(tline(1))
                    p = p+1;
                    pqo.species.name{n,p} = strtrim(tline(1:19));
                    data = sscanf(tline,'%*s %f %*f %*f %*f %*f %*f');
                    pqo.species.molal(n,p) = data(1);
                end
                tline = fgetl(fid);
            end
            
            solids_mass = solids.mass(ismember(solids.str,pqo.solids(n,:)'));
            solids_rho = solids.rho(ismember(solids.str,pqo.solids(n,:)'));
            
             % mass of solid salts (hydrates) in brine
              pqo.sm(n) = sum(pqo.solid_moles(n,~strcmp(pqo.solids(n,:),'Ice(s)'))'.*solids_mass(~strcmp(pqo.solids(n,:),'Ice(s)')),[],'omitnan');
            
            if pqo.sm(n)==0
                pqo.rho_ss(n) = nan;
            else
                ind = find(pqo.solid_moles(n,:)>0 & ~strcmp(pqo.solids(n,:),'Ice(s)'));
                pqo.rho_ss(n) = sum(pqo.solid_moles(n,ind)'.*...
                    solids_mass(ind),[],'omitnan')/sum((pqo.solid_moles(n,ind)'.*...
                    solids_mass(ind))./solids_rho(ind),[],'omitnan');
                if isnan(pqo.rho_ss(n))
                end
            end
            
        end
    end
end
fclose(fid);

if all(isnan(pqo.pressure)) % Pressure only outputs if not 1 atm
    pqo.pressure = ones(size(pqo.pressure));
end

% Derived bulk properties
pqo.br_H2O = pqo.mass_H2O*1000; % Mass of Water in Brine (g)
ind_ice = find(strcmp(pqo.solids(1,:),'Ice(s)'));
pqo.ice_H2O = pqo.solid_moles(:,ind_ice)*solids.mass(strcmp(solids.str,'Ice(s)')); % Mass of Ice (g)

sol_species = pqo.species.name(1,:);

pqo.br_salt = 0;
for p = 1:length(sol_species)
    pqo.br_salt = pqo.br_salt + ...
        pqo.mass_H2O.*(species.mass(strcmp(sol_species{p},species.str))*pqo.species.molal(:,p));
end




% Eutectic
m_crit = 1e-3; % kg
ind_crit =  find(pqo.mass_H2O<m_crit);
if length(ind_crit)==1
    ind_eut = ind_crit;
elseif isempty(ind_crit)
    ind_eut = length(pqo.T);
else
    ind_eut = min(ind_crit);  
end
pqo.eutectic.T = pqo.T(ind_eut);

disp(['Solids precipitating at T = ',num2str(pqo.T(ind_eut)),' K'])
ind = find(pqo.solid_moles(ind_eut,:)>0);
solids_str = pqo.solids(1,ind);
solids_disp = solids_str;
text1 = sprintf('%s\n',solids_disp{:});
fprintf('%s',text1);


end