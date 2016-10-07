% create_Tempest_TileFile_LL2CS.m
% This script converts Tempest LL2CS output to xData and creates a
% tile file that roughly regrids the same as the equivalent GMAO tilefile.

% NOTES:
%  (1) It is intended to create tile files at resolutions that we do not have 
%      from GAMO. 
%  (2) Use compare_Tempest_GMAO_LL2CS.m to see the differences between the
%      Tempest and GMAO tile files for resolutions for which we have both.
%
% Lizzie Lundgren, 10/7/16

clear all; close all

% Tools and work paths
homeDir = '/n/home08/elundgren/GCHP/tools/';
CSGridDir = [homeDir, 'CSGrid/'];
addpath(genpath(CSGridDir));

% Set grid parameters
Nlon = 1440;
Nlat = 720;
NX   = 24;

% Set Tempest LL2CS netcdf info
% NOTE: files created using ~elundgren/GCHP/tools/Tempest/runTempest.sh
%       with unaltered tempestremap copied from /n/home03/zhuangjw/test/.
%       I sourced the GCHP bashrc prior to compiling and running.
tempStr = ['lon' num2str(Nlon) '_lat' num2str(Nlat)];
NXStr = num2str(NX);
tempestPath = [homeDir, 'Tempest/output/' tempStr '_c' NXStr '/'];
tempestFile =  [tempStr '-to-c' NXStr '_MAP.nc'];

% Set Tempest tile file info
lonStr = [repmat('0', 1, 4-length(num2str(Nlon))), num2str(Nlon)];
latStr = [repmat('0', 1, 4-length(num2str(Nlat))), num2str(Nlat)];
cStr   = [repmat('0', 1, 4-length(num2str(NX))),   num2str(NX)];
TileFileName = ['DE' lonStr 'xPE' latStr '_CF' cStr 'x6C.bin']; 
tempestTileFile = [homeDir, 'TileFiles/Tempest/' TileFileName]; 

% Create tile file
xData_temp = readTempest( [tempestPath, tempestFile], Nlon, Nlat, NX);
writeTileFile( tempestTileFile, xData_temp );

% Display contents for validation
xData_new  = displayTileFile( tempestTileFile );

