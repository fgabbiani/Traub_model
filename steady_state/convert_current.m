function [ traub_current_nA ] = convert_current( pr_current )
%convert_current: from Traub to Pinsky-Rinzel model
%   pr_current = convert_current( traub_current )
%
%   where: traub_current = traub model injected somatic current in nA 
%          pr_current = equivalent Pinsky-Rinzel current in uA/cm^2

%fraction of the membrane area attributed to the soma
p = 0.5;

soma_area = 3320; % um^2
basilar_area = 1673;
apical_area = 2188;

n_basilar = 8; 
n_apical = 10;

%total area computed from each compartment area and the number of
%compartments
total_area_um2 = soma_area + n_basilar*basilar_area + n_apical*apical_area;

%conversion factor from um to cm
um2cm = 1e-4; % 1 um = 10^-6 m = 10^-4 cm (1 cm = 10^-2 m)

%total area in cm2
total_area_cm2 = total_area_um2 * um2cm * um2cm;

%conversion factor from nA to uA
nA2uA = 1e-3; % 1 nA = 10^-9 A = 10^-3 uA (1 uA = 10^-3 A)

%scale per unit area and assign fraction attributed to soma as in legend of
%fig. 6 of Pinsky-Rinzel paper
%pr_current = traub_current_uA/(p*total_area_cm2);
traub_current_nA = pr_current*(p*total_area_cm2)/(nA2uA);
end

