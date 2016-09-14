%====================================================
% compare_hp_to_classic.m
%
% Script to do basic comparison between GCHP and GEOS-Chem Classic output
% using S. Eastham's CSGrid matlab toolbox (https://bitbucket.org/gcst/csgrid)
%
% E. Lundgren, 9/6/16
%====================================================
clear all; close all;

global CSGridDir;

% Working directory
workdir = pwd;

% GCHP home directory
homedir = '/n/regal/jacob_lab/elundgren/GCHP/';

% GC run directories (GCHP and Classic)
dir_hp = [homedir, 'testruns/gchp_4x5_tropchem/OutputDir/'];
dir_cc = [homedir, 'testruns/geosfp_4x5_tropchem_classic/'];  

% Output netcdf files to compare
file_hp_cs = [dir_hp, 'GCHP.center.20130701.nc4'];
file_hp_ll = [dir_hp, 'GCHP.regrid.20130701.nc4'];
file_cc    = [dir_cc, 'GEOSCHEM_Diagnostics_Hrly.201307010100.nc'];

% Add necessary paths
CSGridDir = '/n/regal/jacob_lab/elundgren/GCHP/tools/CSGrid';
addpath(genpath(CSGridDir));

% Define species of interest for each file
spclist = {'O3','NO','CO','NO2'};

% Define level of interest, where index increases starting at the surface
% NOTE: This requires GCHP data to be vertically inverted which is done below
lev = 1;

% Declare structures to store data
hp_cs = struct;
hp_ll = struct;
cc    = struct;
spc   = struct;

% Read raw data for species of interest
for i = 1:length(spclist);

  % Make local variable for species
  var = spclist{i};

  % Declare structures
  hp_cs.(var) = struct;
  hp_ll.(var) = struct;
  cc.(var) = struct;

  % Read raw data for that variable
  hp_cs.(var).data = ncread(file_hp_cs, ['TRC_', var]);
  hp_ll.(var).data = ncread(file_hp_ll, ['TRC_', var]);
  cc.(var).data    = ncread(file_cc,    ['TRACER_CONC_', var]);
  
  % ***IMPORTANT!!!*** Vertically invert the data (GCHP only)
  hp_cs.(var).data = flip(hp_cs.(var).data,3);
  hp_ll.(var).data = flip(hp_ll.(var).data,3);

  % Store min and max for each dataset
  hp_cs.(var).min = min(hp_cs.(var).data(:));
  hp_cs.(var).max = max(hp_cs.(var).data(:));
  hp_ll.(var).min = min(hp_ll.(var).data(:));
  hp_ll.(var).max = max(hp_ll.(var).data(:));
  cc.(var).min    = min(cc.(var).data(:));
  cc.(var).max    = max(cc.(var).data(:));

  % store min and max for each species
  spc.(var) = struct;
  spc.(var).min = min([hp_cs.(var).min   ...
                       hp_ll.(var).min   ...
                       cc.(var).min ]);
  spc.(var).max = max([hp_cs.(var).max   ... 
		       hp_ll.(var).max   ...
		       cc.(var).max ]);
end

%-----------------------
% Plot surface ozone
%-----------------------

var = 'O3';

% Plot GCHP lat/lon O3 surface data
figure(1)
plotGrid(hp_ll.(var).data(:,:,1),'layer','projection','miller');
colorbar;
%caxis([spc.O3.min spc.O3.max]);
title(['GCHP 4x5 regridded, surface ', var]);


% Plot GEOS-Chem Classic O3 surface data
figure(2)
plotGrid(cc.(var).data(:,:,1),'layer','projection','miller');
colorbar;
%caxis([spc.O3.min spc.O3.max]);
title(['GEOS-Chem 4x5, surface ', var]);

%% GCHP cubed-sphere
%figure(3)
%titletext = 'GCHP cubed-sphere';
%plotCSLayer(hp_cs.(var).data(:,:,1),'projection','miller');
%colorbar;
%%caxis([spc.O3.min spc.O3.max]);




