% compare_Tempest_GMAO_LL2CS.m
% This script converts Tempest LL2CS output to xData and tests whether it
% regrids lat/lon to cubed sphere roughly the same as a GMAO tilefile
% for the same resolutions. The comparison is made by the user by visually
% inspecting output figures.

clear all; close all

%=======================================================================
% Test case setup (edit this section to match your local settings/files)
%=======================================================================
% Tools and work paths
homeDir = '/n/home08/elundgren/GCHP/tools/';
CSGridDir = [homeDir, 'CSGrid/'];

% Choose what test to run (options are defined in gridParams)
testid = 5;

% Set test case grid parameters
gridParams = [[72  46  48 ]; % test 1 lon, lat, cs
              [360 180 24 ]; % test 2 lon, lat, cs
              [360 180 48 ]; % etc
              [360 180 90 ];
              [360 180 180]];

% Set grid parameters for the chosen test case
Nlon = gridParams( testid, 1 );
Nlat = gridParams( testid, 2 );
NX   = gridParams( testid, 3 );

% Set Tempest LL2CS netcdf info
% NOTE: files created using ~elundgren/GCHP/tools/Tempest/runTempest.sh
%       with unaltered tempestremap copied from /n/home03/zhuangjw/test/.
%       I sourced the GCHP bashrc prior to compiling and running.
tempStr = ['lon' num2str(Nlon) '_lat' num2str(Nlat)];
NXStr = num2str(NX);
tempestPath = [homeDir, 'Tempest/output/' tempStr '_c' NXStr '/'];
tempestFile =  [tempStr '-to-c' NXStr '_MAP.nc'];

% Set GMAO tilefile info
% NOTE: files copied from /n/regal/jacob_lab/mslong/ForSeb/
gmaoPath = [homeDir, 'TileFiles/GMAO/'];
lonStr = [repmat('0', 1, 4-length(num2str(Nlon))), num2str(Nlon)];
latStr = [repmat('0', 1, 4-length(num2str(Nlat))), num2str(Nlat)];
cStr   = [repmat('0', 1, 4-length(num2str(NX))),   num2str(NX)];
gmaoFile = ['DE' lonStr 'xPE' latStr '_CF' cStr 'x6C.bin']; 

% Define Tempest tile file to write (you will be asked if you want to save)
% NOTE: use same name as GMAO tilefile for ready compatibility with ExtData
tempestTileFile = [homeDir, 'TileFiles/Tempest/' gmaoFile]; 

addpath(genpath(CSGridDir));

% To save plots:
plotsdir = '~elundgren/GCHP/tools/Plots/Tempest_vs_GMAO_tilefiles/';
addpath(genpath('~elundgren/matlab/export_fig'));
addpath(genpath('~elundgren/matlab/b2r')); % for blue to red colormap

%======================================================
%  Generate fake reference CS data for testing
%======================================================
CSdata_fake = zeros(NX,NX*6);

% a "F"-like graph to check the cube-face orientation. Seems silly.. 
F = zeros(NX,NX);
Fwidth = fix(NX/15);
minX = fix(NX/3);
maxX = fix(NX*2/3);
minY = fix(NX/4);
midY = fix(NX/2);
maxY = fix(NX*3/4);
F( minX : minX + 2*Fwidth, minY : maxY                   ) = 1;
F( minX : maxX,            midY - Fwidth : midY + Fwidth ) = 1;
F( minX : maxX,                   maxY - Fwidth : maxY + Fwidth ) = 1;

% Set background value to the panel number, include F in 'fake' data
for i=1:6
    CSdata_bg(:, i*NX-NX+1 : i*NX) = i * ones(NX,NX); 
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

% Save figure
filename = ['CSorig_', tempStr, '_c', NXStr]; 
eval(['export_fig ', [plotsdir, filename], '.png -nocrop']); 

%======================================================
% Use existing GMAO tilefile for regridding
%======================================================
xData_gmao = readTileFile([ gmaoPath, gmaoFile]);

% regrid cubed sphere to lat/lon using GMAO tile file
[ LLdata_ref ] = regridAndPlot( CSdata_fake , xData_gmao, 'CS2LL', ...
			        'original cs --gmao--> lat/lon'); 

% regrid lat/lon back to cubed sphere (reference)
[ CSdata_GMAO ] = regridAndPlot( LLdata_ref, xData_gmao, 'LL2CS', ...
			         'lat/lon --gmao--> cs (ref)');

% Save figure
filename = ['CSgmao_', tempStr, '_c', NXStr]; 
eval(['export_fig ', [plotsdir, filename], '.png -nocrop']); 

%==================================================
% Read/transform Tempest LL2CS netcdf info to xData
%==================================================
xData_temp = readTempest( [tempestPath, tempestFile], Nlon, Nlat, NX);

%==================================================
% Regrid lat/lon to CS to compare with reference
%==================================================
[ CSdata_tempest ] = regridAndPlot( LLdata_ref, xData_temp, 'LL2CS', ...
				   'lat/lon --tempest--> cs (after fixes)');

% Save figure
filename = ['CStempest_', tempStr, '_c', NXStr]; 
eval(['export_fig ', [plotsdir, filename], '.png -nocrop']); 

%=======================================================================
% Plot the regrid differences (Tempest - GMAO), normalized by background
%=======================================================================
CSdata_diff = ( CSdata_tempest - CSdata_GMAO ) ./ CSdata_bg;
[ CSdata_diff ] = regridAndPlot( CSdata_diff, xData_temp, 'LL2CS', ...
				 [num2str(Nlon) 'x' num2str(Nlat), ...
				  ' -> ', NXStr, ...
				  ': ( Tempest - GMAO ) / background'], ...
				 false, true); 
		                 % Last two args are optional:
                                 % passing 1st arg as false skips regridding,
		                 % passing 2nd arg as true uses blue to red cmap

% Save figure
filename = ['GMAO_vs_Tempest_', tempStr, '_c', NXStr]; 
eval(['export_fig ', [plotsdir, filename], '.png -nocrop']); 

%====================================
% Create tilefile with Tempest xData?
%====================================
prompt = 'Do you want to create a new tile file from the Tempest data? (y/n)\n';
answer = input(prompt,'s');
if strcmp(answer,'y');
   writeTileFile( tempestTileFile, xData_temp );

   %==========================================================
   % Display GMAO and new Tempest tilefile info for validation
   % NOTE: Tempest and GMAO numPoints differ in xData.
   %       This is to be expected. However, they should be
   %       approximately the same.
   %==========================================================
   xData_gmao = displayTileFile( [ gmaoPath, gmaoFile] );
   xData_new  = displayTileFile( tempestTileFile );
elseif strcmp(answer,'n');
   fprintf('No tile file created\n');		 
else
   error('not a valid answer');
end









