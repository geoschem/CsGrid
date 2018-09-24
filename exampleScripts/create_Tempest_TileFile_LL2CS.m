% create_Tempest_TileFile_LL2CS.m
% This script converts Tempest LL2CS nc output to binary for read by MAPL
%
% ~ Lizzie Lundgren

clear all; close all

% Tools and work paths
CSGridDir = '/n/home08/elundgren/GC/matlab/CSGrid/';
addpath(genpath(CSGridDir));

% Set grid parameters
Nlon = 1250; % longitude size
Nlat = 700;  % latitude size
NX   = 24;   % cubed sphere res
dateline = 'UU'; % DE or DC for global; UU for regional
polar = 'UU'; % PE or PC for global; UU for regional

% Get strings for filenames
lonStr = [repmat('0', 1, 4-length(num2str(Nlon))), num2str(Nlon)];
latStr = [repmat('0', 1, 4-length(num2str(Nlat))), num2str(Nlat)];
cStr   = [repmat('0', 1, 4-length(num2str(NX))),   num2str(NX)];

% Define tile name prefix
tfPrefix = [dateline lonStr 'x' num2str(polar) latStr '_CF' cStr 'x6C']

% Define Tempest netcdf file path (to be read)
ncFileName =  [tfPrefix '.nc']
ncFileDir = '~elundgren/GC/regridding/tempestremap/TileFiles/';
ncFilePath = [ncFileDir ncFileName]

% Define Tempest binary tile file path (to be written)
binFileName = [tfPrefix '.bin']
binFileDir = '~elundgren/GC/regridding/tempestremap/TileFiles/'; 
binFilePath = [binFileDir binFileName] 

% Read netcdf tile file
xData_nc = readTempest( ncFilePath );

% Create binary tile file
writeTileFile( binFilePath, xData_nc );

% Display contents for validation
xData_bin  = displayTileFile( binFilePath );

