% testTempest_LL2CS.m
% This script converts Tempest LL2CS output to xData and tests whether it
% regrids lat/lon to cubed sphere roughly the same as a GMAO tilefile
% for the same resolutions.

clear all; close all

%==========================================
%  Configurable parameters
%==========================================

% Set grid parameters
NX=48;
Nlon=72;
Nlat=46;

% Tempest LL2CS netcdf info (make sure file matches grid params!)
tempestPath = '/n/home08/elundgren/GCHP/tools/tilefile_from_Jaiwei/tilefile/';
tempestFile =  '4x5-to-c48_MAP.nc';

% GMAO tilefile info (make sure file matches grid params!)
gmaoPath = tempestPath;
gmaoFile =  'DE0072xPE0046_CF0048x6C.bin';

% Add matlab paths
CSGridDir = '/n/home08/elundgren/GCHP/tools/CSGrid/';
addpath(genpath(CSGridDir));

%======================================================
%  Generate reference CS data for testing (do not edit)
%======================================================
CSdata_fake = zeros(NX,NX*6);

% a "F"-like graph to check the cube-face orientation. Seems silly.. 
F = zeros(NX,NX);
Fwidth = 2; 
minX = fix(NX/3);
maxX = fix(NX*2/3);
minY = fix(NX/4);
midY = fix(NX/2);
maxY = fix(NX*3/4);
F( minX - Fwidth : minX + Fwidth, minY : maxY                   ) = 1;
F( minX : maxX,                   midY - Fwidth : midY + Fwidth ) = 1;
F( minX : maxX,                   maxY - Fwidth : maxY + Fwidth ) = 1;

% Set background value to the panel number
for i=1:6
    CSdata_fake(:, i*NX-NX+1 : i*NX) = i * ( ones(NX,NX) - F ); 
end

% plot
figure;
for i=1:6
    subplot(2,3,i)
    surf(CSdata_fake(:,i*NX-NX+1:i*NX)','EdgeColor','None');
    colorbar
    caxis([0,6])
    xlim([1,NX])
    ylim([1,NX])
    view(2); 
    title(sprintf('panel %d',i));
end   
suptitle('original cs')

% Read GMAO tilefile
xData_gmao = readTileFile([ gmaoPath, gmaoFile]);

%  regrid cubed sphere to lat/lon using GMAO tile file
[ LLdata_ref ] = regridAndPlot( CSdata_fake , xData_gmao, 'CS2LL', ...
			       'original cs --gmao--> lat/lon');

%  regrid lat/lon back to cubed sphere (reference)
[ CSdata_GMAO ] = regridAndPlot( LLdata_ref, xData_gmao, 'LL2CS', ...
			       'lat/lon --gmao--> cs (ref)');

%==================================================
% Read/transform Tempest LL2CS netcdf info to xData
%==================================================

% Get tempest LL2CS netcdf file information into xData
xData_temp = getTempest_LL2CS_xData( [tempestPath, tempestFile], ...
				      xData_gmao(1).Name, ...
				      xData_gmao(2).Name, Nlat, Nlon, NX);

% NOTE: eventually add all the below fixes to the function called above.

% Regrid lat/lon to CS using raw Tempest xData
[ CSdata_tempest ] = regridAndPlot( LLdata_ref, xData_temp, 'LL2CS', ...
				   'lat/lon --tempest--> cs (raw)');	

% Correct longitude shift
% NOTE: need to make this offset automatic from resolution
xData_temp(1).II = shiftIndexes( xData_temp(1).II, -2, Nlon );

% Swap faces
xData_temp(2) = swapCSFaces( xData_temp(2), [4 5 1 2 6 3]);

% Flip face 6 in both dimensions
xData_temp(2) = flipCSFaces( xData_temp(2), 6, '2d');

% Transpose and flip (II only) faces 3-5
xData_temp(2) = transposeCSFaces( xData_temp(2), [3 4 5]);
xData_temp(2) = flipCSFaces( xData_temp(2), [3 4 5], 'II');

%==================================================
% Regrid lat/lon to CS to check if matches ref
%==================================================
[ CSdata_tempest ] = regridAndPlot( LLdata_ref, xData_temp, 'LL2CS', ...
				   'lat/lon --tempest--> cs (after fixes)');







