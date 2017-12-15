cd% create_Tempest_TileFile_LL2CS.m
% This script converts Tempest LL2CS nc output to binary for read by MAPL
%
% ~ Lizzie Lundgren

clear all; close all

% Tools and work paths
CSGridDir = '/n/home08/elundgren/GC/matlab/CSGrid/';       % CSGrid full path
TempestFileDir = CSGridDir;  % Directory where your tempest output is
TileFileDir = CSGridDir;     % Where you want the binary tilefile output
addpath(genpath(CSGridDir));

% Set grid parameters
Nlon = 360;
Nlat = 180;
NX   = 24;
dateline = 'DC';
polar = 'PC';

% Set Tempest LL2CS netcdf info
cellType = [dateline '_' polar];
tempStr = ['lon' num2str(Nlon) '_lat' num2str(Nlat)];
NXStr = num2str(NX);
TempestFile =  [tempStr '-to-c' NXStr '_MAP_' cellType '.nc'];

% Set Tempest tile file info
lonStr = [repmat('0', 1, 4-length(num2str(Nlon))), num2str(Nlon)];
latStr = [repmat('0', 1, 4-length(num2str(Nlat))), num2str(Nlat)];
cStr   = [repmat('0', 1, 4-length(num2str(NX))),   num2str(NX)];
TileFileName = [dateline lonStr 'x' polar latStr '_CF' cStr 'x6C.bin']; 
TempestTileFile = [TileFileDir TileFileName]; 

% Create tile file
xData_temp = readTempest( [TempestFileDir, TempestFile] );
writeTileFile( TempestTileFile, xData_temp );

% Display contents for validation
xData_new  = displayTileFile( TempestTileFile );

