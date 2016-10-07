% replicateTileFile.m
%   Script illustrating how to replicate an existing tilefile by reading
%   the binary, storing the data, writing to a new binary file,
%   and reading the new file for validation. This script can be adapted
%   to create a new tilefile from a new xData array structure.
% 
% Lizzie Lundgren, 9/26/16 

clear all 
close all

% Add paths
CSGridDir = '/n/regal/jacob_lab/elundgren/GCHP/tools/CSGrid';
addpath(genpath(CSGridDir));

% Define existing tilefile (to read)
tf_path_orig = [CSGridDir, '/GridData/TileFiles/'];
tf_name_orig = 'DC0144xPC0091_CF0024x6C.bin';
tf_file_orig = [tf_path_orig, tf_name_orig];

% Define new tilefile (to write)
tf_path_new = '~elundgren/GCHP/tools/';
tf_name_new = 'DC0144xPC0091_CF0024x6C_new.bin';
tf_file_new = [tf_path_new, tf_name_new];

% Read original tilefile, and store data in xData_orig
xData_orig = displayTileFile(tf_file_orig);

% Create new tilefile that replicates the old by writing xData_orig
s = writeTileFile( tf_file_new, xData_orig ); 
assert ( s~=1, 'Error in writeTileFile');
if s == 0;
  fprintf('Success!\n');

  % Read new tilefile for validation
  xData_new = displayTileFile(tf_file_new);

elseif s == 2;
  fprintf('Exiting replicateTileFile since tilefile already exists.\n');
end  