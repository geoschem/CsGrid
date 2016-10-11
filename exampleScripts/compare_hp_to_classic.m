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

%-------------------
% Configurables
%-------------------

% Define species of interest (may include O3, NO, CO, NO2)
spclist = {'O3'};

% Define level of interest, where index increases starting at the surface.
% This must be scalar. NOTE: This requires GCHP data to be vertically 
% inverted which is done below.
lev = 1;

% Define run directory prefixes for runs you did. Use the same prefix
% all GCHP runs, and the same prefix for GC classic runs.
hp_run_prefix = 'hp_4x5_trop_';
cc_run_prefix = 'classic_4x5_trop_';

% Define run directory suffixes for runs you did (same for GCHP and classic)
% These are also used as structure field names so choose wisely!
% For example, this script is set for the following runs:
%   std    : out of the box settings
%   noET   : no emissions or transport
%   noETDD : no emissions, transport, or drydep
% This corresponds to GCHP run directories hp_4x5_trop_std, 
% hp_4x5_trop_noET, etc, and GC classic run directories classic_4x5_trop_std,
% classic_4x5_trop_noET, etc. 
run_suffix = {'DDonly', 'allOff'};

% Choose which model outputs to include. This is where you can isolate
% comparing only GCHP, or only classic.
gchp_cubedsphere_on = true;
gchp_latlon_on      = true;
classic_on          = true;

% Choose whether to plot data (will plot all runs, all species, at level lev)
% NOTE: this can be slow. Best compare without plots first to look
% at metrics.
plots_on = true;

% Choose what metrics to compare between runs. Currently min, max, or mean.
metric = 'mean';

% Testruns directory, where all your run directories are stored
testdir = '/n/regal/jacob_lab/elundgren/GCHP/testruns/';

% Your local CSGrid repository location
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
% Loop over runs and species to read in data
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
      hp_cs.(run).(spc).min  = min(min(hp_cs.(run).(spc).data(:,:,lev)));
      hp_cs.(run).(spc).max  = max(max(hp_cs.(run).(spc).data(:,:,lev)));
      hp_cs.(run).(spc).mean = mean(mean(hp_cs.(run).(spc).data(:,:,lev)));
    end  

    % GCHP runs, lat/lon output
    if gchp_latlon_on;
      hp_ll.(run).(spc)      = struct;
      hp_ll.(run).(spc).data = ncread(fpath_hp_ll{i}, ['TRC_', spc]);
      hp_ll.(run).(spc).data = flip(hp_ll.(run).(spc).data,3);      
      hp_ll.(run).(spc).min  = min(min(hp_ll.(run).(spc).data(:,:,lev)));
      hp_ll.(run).(spc).max  = max(max(hp_ll.(run).(spc).data(:,:,lev)));
      hp_ll.(run).(spc).mean = mean(mean(hp_ll.(run).(spc).data(:,:,lev)));
    end  

    % GEOS-Chem classic runs
    if classic_on;
      cc.(run).(spc)      = struct;      
      cc.(run).(spc).data = ncread(fpath_cc{i}, ['TRACER_CONC_', spc]);
      cc.(run).(spc).min  = min(min(cc.(run).(spc).data(:,:,lev)));
      cc.(run).(spc).max  = max(max(cc.(run).(spc).data(:,:,lev)));
      cc.(run).(spc).mean = mean(mean(cc.(run).(spc).data(:,:,lev)));
    end  
  end
end

%---------------------------------------------------------
% Get min and max for all runs per species for comparison
%---------------------------------------------------------
cminval = zeros(1,length(spclist));
cmaxval = zeros(1,length(spclist));
for i =1:length(run_suffix);
  run = run_suffix{i};
  for j = 1:length(spclist);
    spc = spclist{j};
    cminval(j)=1e30;
    cmaxval(j)=1e-30;

    if gchp_cubedsphere_on;
      if hp_cs.(run).(spc).min < cminval(j);
	cminval(j) = hp_cs.(run).(spc).min;
      end
      if hp_cs.(run).(spc).max > cmaxval(j);
	cmaxval(j) = hp_cs.(run).(spc).max;
      end
    end  
  
    % GCHP runs, lat/lon output
    if gchp_latlon_on;
      if hp_ll.(run).(spc).min < cminval(j);
        cminval(j) = hp_ll.(run).(spc).min;
      end
      if hp_ll.(run).(spc).max > cmaxval(j);
        cmaxval(j) = hp_ll.(run).(spc).max;
      end
    end  
  
    % classic runs
    if classic_on;
      if cc.(run).(spc).min < cminval(j);
        cminval(j) = cc.(run).(spc).min;
      end
      if cc.(run).(spc).max > cmaxval(j);
        cmaxval(j) = cc.(run).(spc).max;
      end
    end 
  end
end

%------------------------------------------------------
% Loop over runs and species and plot
%------------------------------------------------------
if plots_on;
  for i =1:length(run_suffix);
    run = run_suffix{i};
    for j = 1:length(spclist);
      spc = spclist{j};
  
      % GCHP runs, cubed sphere output
      if gchp_cubedsphere_on;
        figure(h)
        titletext = 'GCHP cubed-sphere';
        plotCSLayer(hp_cs.(run).(spc).data(:,:,lev),'projection','miller');
        colorbar;
        caxis([cminval(j) cmaxval(j)]);
        % NOTE: title not currently possible plotCSLayer?
        h = h + 1;
      end  
  
      % GCHP runs, lat/lon output
      if gchp_latlon_on;
        figure(h)
        plotGrid(hp_ll.(run).(spc).data(:,:,lev), ...
        	       'layer','projection','miller');
        colorbar;
        caxis([cminval(j) cmaxval(j)]);
        title(['GCHP 4x5 regridded, GC level ', num2str(lev), ...
        	     ', ', spc, ', ', run]);
        h = h + 1;
      end  
  
      % classic runs
      if classic_on;
        figure(h)
        plotGrid(cc.(run).(spc).data(:,:,lev),'layer','projection','miller');
        colorbar;
        caxis([cminval(j) cmaxval(j)]);
        title(['Classic 4x5, GC level ', num2str(lev), ...
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
  fprintf('\n%s values for %s, level %d\n',metric,spc,lev);
  for j =1:length(run_suffix);
    run = run_suffix{j};
    fprintf(' Run %s\n',run); 
    if gchp_cubedsphere_on; 
      fprintf('   cubed sphere: %12.8d\n',hp_cs.(run).(spc).(metric)); 
    end
    if gchp_latlon_on;      
      fprintf('   cs regridded: %12.8d\n',hp_ll.(run).(spc).(metric)); 
    end
    if classic_on;          
      fprintf('   classic:      %12.8d\n',cc.(run).(spc).(metric)); 
    end
  end
end

%---------------------------------------------------------------- 
% Get difference of HP (regridded) and classic per spc/run combo
%----------------------------------------------------------------
if gchp_latlon_on && classic_on;
  for i =1:length(run_suffix);
    run = run_suffix{i};
    for j = 1:length(spclist);
      spc = spclist{j};

      % calculate data difference and ratio
      datadiff = hp_ll.(run).(spc).data(:,:,lev) ...
                 - cc.(run).(spc).data(:,:,lev);
      dataratio = hp_ll.(run).(spc).data(:,:,lev) ...
                  ./cc.(run).(spc).data(:,:,lev);
  
      % plot the difference
      if plots_on;
        figure(h)
        plotGrid(datadiff,'layer','projection','miller');
        title([run ': HP - classic, level ', ...
	       num2str(lev), ', ', spc]);
        colorbar;
        h = h + 1;
      end  

      % plot the ratio
      if plots_on;
        figure(h)
        plotGrid(dataratio,'layer','projection','miller');
        title([run ': HP / classic, level ', ...
	       num2str(lev), ', ', spc]);
        colorbar;
        h = h + 1;
      end  
    end
  end
end

%----------------------------------------------------
% Get HP differences of run 1 vs run2 per species
%----------------------------------------------------
run1 = run_suffix{1};
run2 = run_suffix{2};
for j = 1:length(spclist);
  spc = spclist{j};

  % gchp lat/lon (regridded0
  if gchp_latlon_on;
     datadiff = hp_ll.(run1).(spc).data(:,:,lev) ...
                - hp_ll.(run2).(spc).data(:,:,lev);
     dataratio = hp_ll.(run1).(spc).data(:,:,lev) ...
                 ./hp_ll.(run2).(spc).data(:,:,lev);
     if plots_on;
       figure(h)
       plotGrid(datadiff,'layer','projection','miller');
       title(['HP: ', run1, ' - ', run2, ...
     	     ', level ', num2str(lev), ', ', spc]);
       colorbar;
       h = h + 1;

       figure(h)
       plotGrid(dataratio,'layer','projection','miller');
       title(['HP: ', run1, ' / ', run2, ...
     	     ', level ', num2str(lev), ', ', spc]);
       colorbar;
       h = h + 1;
     end
  end  

  % classic
  if classic_on;
     datadiff = cc.(run1).(spc).data(:,:,lev) ...
                - cc.(run2).(spc).data(:,:,lev);
     dataratio = cc.(run1).(spc).data(:,:,lev) ...
                 ./cc.(run2).(spc).data(:,:,lev);
     if plots_on;
       figure(h)
       plotGrid(datadiff,'layer','projection','miller');
       title(['Classic: ', run1, ' - ', run2, ...
     	     ', level ', num2str(lev), ', ', spc]);
       colorbar;
       h = h + 1;

       figure(h)
       plotGrid(dataratio,'layer','projection','miller');
       title(['Classic: ', run1, ' / ', run2, ...
     	     ', level ', num2str(lev), ', ', spc]);
       colorbar;
       h = h + 1;
     end 
  end
end





