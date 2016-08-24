function [ kgGrid ] = vvToMass( vvGrid, pEdge_hPa, gridArea_m2, kgPerMol )
%VVTOMASS Convert v/v to kg
%   vvGrid:         Mixing ratios (mol X/mol air)
%   pEdge_hPa:      3D array of box edge pressures (hPa)
%   gridArea_m2:    2D array of grid box areas (m2)
%   kgPerMol:       Molar mass of X (kg/mol)

persistent molMassAir;
if isempty(molMassAir)
    load constData MWAir_kg;
    molMassAir = MWAir_kg;
end

airMass = calcAirMass(pEdge_hPa,gridArea_m2);
massMR = vvGrid.*(kgPerMol./molMassAir);
kgGrid = massMR.*airMass;

end

