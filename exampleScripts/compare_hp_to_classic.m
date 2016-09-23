%====================================================
% compare_hp_to_classic.m
%
% Script to do basic comparison between GCHP and GEOS-Chem Classic outputs
% using S. Eastham's CSGrid matlab toolbox (https://bitbucket.org/gcst/csgrid)
%
% E. Lundgren, 9/6/16
%====================================================
clear all; close all;

global CSGridDir
format long

%-------------------
% Configurables
%-------------------

% Define species of interest (may include O3, NO, CO, NO2)
spclist = {'O3', 'NO'};

% Define level of interest, where index increases starting at the surface.
% This must be scalar. NOTE: This requires GCHP data to be vertically 
% inverted which is done below.
lev = 1;

% Define run directory suffixes for runs you did (same for GCHP and classic)
% These are also used as structure field names so choose wisely!
% My current runs include:
%   std    : out of the box settings
%   noET   : no emissions or transport
%   noETDD : no emissions, transport, or drydep
run_suffix = {'noETDD', 'std', 'noET', 'noETDDC'};

% Choose which ones to include
gchp_cubedsphere_on = true;
gchp_latlon_on      = true;
classic_on          = true;

% Choose whether to plot data (will plot all runs, all species, at level lev)
plots_on = false;

% Define run directory prefixes for runs you did. Use the same prefix
% all GCHP runs, and the same prefix for GC classic runs.
hp_run_prefix = 'hp_4x5_trop_';
cc_run_prefix = 'classic_4x5_trop_';

% Testruns directory, where all run directories are stored
testdir = '/n/regal/jacob_lab/elundgren/GCHP/testruns/';

% Local CSGrid repository location
CSGridDir = '/n/regal/jacob_lab/elundgren/GCHP/tools/CSGrid';

%---------------------------------------------------------
% The rest... don't need to change unless adding features
%---------------------------------------------------------

% Add CSGrid path
addpath(genpath(CSGridDir));

% Output filenames
fn_hp_cs = 'GCHP.center.20130701.nc4';
fn_hp_ll = 'GCHP.regrid.20130701.nc4';
fn_cc    = 'GEOSCHEM_Diagnostics_Hrly.201307010100.nc';

% Define output file paths
fpath_hp_cs = cell(3,1); % cubed-sphere GCHP output
fpath_hp_ll = cell(3,1); % lat/lon GCHP output
fpath_cc    = cell(3,1); % lat/lon GEOS-Chem classic output
for i = 1:length(run_suffix);
  fpath_hp_cs{i} = [testdir hp_run_prefix run_suffix{i} '/OutputDir/' fn_hp_cs];
  fpath_hp_ll{i} = [testdir hp_run_prefix run_suffix{i} '/OutputDir/' fn_hp_ll];
  fpath_cc{i}    = [testdir cc_run_prefix run_suffix{i} '/' fn_cc];
end

% Declare structures to store data
if gchp_cubedsphere_on; hp_cs = struct; end % cubed sphere GCHP output
if gchp_latlon_on;      hp_ll = struct; end % lat/lon GCHP output (regridded)
if classic_on;          cc    = struct; end % GC classic output

% initialize figure number
h = 1;

%------------------------------------------------------
% Loop over runs and species to read/set data and plot
%------------------------------------------------------
for i =1:length(run_suffix);
  
  % Make local variable for this run configuration
  run = run_suffix{i};

  % Declare structures (fast so do for all)
  if gchp_cubedsphere_on; hp_cs.(run) = struct; end
  if gchp_latlon_on;      hp_ll.(run) = struct; end
  if classic_on;          cc.(run)    = struct; end
 
  % Read raw data for species of interest
  for j = 1:length(spclist);
  
    % Make local variable for species
    spc = spclist{j};

    % GCHP runs, cubed sphere output
    if gchp_cubedsphere_on;
      hp_cs.(run).(spc)      = struct;
      hp_cs.(run).(spc).data = ncread(fpath_hp_cs{i}, ['TRC_', spc]);
      hp_cs.(run).(spc).data = flip(hp_cs.(run).(spc).data,3);
      hp_cs.(run).(spc).min  = min(hp_cs.(run).(spc).data(:));
      hp_cs.(run).(spc).max  = max(hp_cs.(run).(spc).data(:));
      hp_cs.(run).(spc).mean = mean(hp_cs.(run).(spc).data(:));
      if plots_on;
         figure(h)
         titletext = 'GCHP cubed-sphere';
         plotCSLayer(hp_cs.(run).(spc).data(:,:,lev),'projection','miller');
         colorbar;
         %caxis([spc.O3.min spc.O3.max]);
         % NOTE: title not currently possible plotCSLayer?
         h = h + 1;
      end
    end  

    % GCHP runs, lat/lon output
    if gchp_latlon_on;
      hp_ll.(run).(spc)      = struct;
      hp_ll.(run).(spc).data = ncread(fpath_hp_ll{i}, ['TRC_', spc]);
      hp_ll.(run).(spc).data = flip(hp_ll.(run).(spc).data,3);      
      hp_ll.(run).(spc).min  = min(hp_ll.(run).(spc).data(:));
      hp_ll.(run).(spc).max  = max(hp_ll.(run).(spc).data(:));
      hp_ll.(run).(spc).mean = mean(hp_ll.(run).(spc).data(:));
      if plots_on;
        figure(h)
        plotGrid(hp_ll.(run).(spc).data(:,:,lev), ...
        	       'layer','projection','miller');
        colorbar;
        %caxis([spc.O3.min spc.O3.max]);
        title(['GCHP 4x5 regridded, GC level ', num2str(lev), ...
        	     ', ', spc, ', ', run]);
        h = h + 1;
      end
    end  

    % GEOS-Chem classic runs
    if classic_on;
      cc.(run).(spc)      = struct;      
      cc.(run).(spc).data = ncread(fpath_cc{i}, ['TRACER_CONC_', spc]);
      cc.(run).(spc).min  = min(cc.(run).(spc).data(:));
      cc.(run).(spc).max  = max(cc.(run).(spc).data(:));
      cc.(run).(spc).mean = mean(cc.(run).(spc).data(:));
      if plots_on;
        figure(h)
        plotGrid(cc.(run).(spc).data(:,:,lev),'layer','projection','miller');
        colorbar;
        %caxis([spc.O3.min spc.O3.max]);
        title(['GEOS-Chem 4x5, GC level ', num2str(lev), ...
        	     ', ', spc, ', ', run]);
        h = h + 1;
      end
    end  
  end
end

%-------------------------------------- 
% Print out mean values for comparison
%--------------------------------------
for i = 1:length(spclist);
  spc = spclist{i};
  fprintf('\nMean values for %s, level %d\n',spc,lev);
  for j =1:length(run_suffix);
    run = run_suffix{j};
    fprintf(' Run %s\n',run); 
    if gchp_cubedsphere_on; 
      fprintf('   cubed sphere: %d\n',hp_cs.(run).(spc).mean); 
    end
    if gchp_latlon_on;      
      fprintf('   cs regridded: %d\n',hp_ll.(run).(spc).mean); 
    end
    if classic_on;          
      fprintf('   classic:      %d\n',cc.(run).(spc).mean); 
    end
  end
end





