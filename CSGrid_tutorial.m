%====================================================
% gchp_analysis.m
%   Script to investigate gchp output and compare
%   with GEOS-Chem Classic
%====================================================

clear all; close all;

global CSGridDir

WORKDIR = pwd;
CSGridDir = '/n/regal/jacob_lab/elundgren/GCHP/tools/CSGrid';
addpath(genpath(CSGridDir));

%-----------------------
% Seb's example code
%-----------------------

%----------------------
% Look at cubed sphere
%----------------------

% Set a dummy variable of size c24 cubed sphere grid
% to a constant value for each face and layer
testData = zeros(24,144,72);
for iFace=1:6
   for iLayer = 1:72
      testData(:,(1:24)+(iFace-1)*24,iLayer) = iFace*iLayer;
   end
end

% Plot variable in native format (c24)
% NOTE: I modified plotCSLayer.m slightly to make the
% the file path correct relative to current directory
figure(1)
plotCSLayer(testData(:,:,1),'projection','flat');
title('6 Faces of the Cubed Sphere Grid')

%----------------
% Look at 2x25 
%----------------

% Regrid from c24 to 2x25
LLData2 = regridConservative(testData,genGridSpec('gmao2x25','geos5'));

% Plot the regridded variable
figure(2)
plotGrid(LLData2,'zonal');
title('Data regridded from c24 to 2x25, Zonal')

% Extract one layer of the plot
figure(3)
plotGrid(LLData2(:,:,1),'layer');
title('Data regridded from c24 to 2x25, Level 1')

% Compare lat-lon vs CS framework grid areas (compare sums!)
gSpec = genGridSpec('gmao2x25','geos5'); 
% now have gSpec.latEdge, gSpec.lonEdge, gSpec.gridArea
ll_total_area = sum(gSpec.gridArea(:))
[lonEdge,latEdge] = calcCSGrid(24); % cubed-sphere
gridArea = calcCSArea(lonEdge,latEdge); % cubed-sphere
cs_total_area = sum(gridArea(:))
pct_diff_cs_vs_2x25 = (ll_total_area - cs_total_area) / ll_total_area * 100

%%----------------
%% Look at 4x5 
%%----------------
%
% This code currently does not work due to missing tile file
%
%% Regrid from c24 to 4x5
%LLData4 = regridConservative(testData,genGridSpec('gmao4x5','geos5'));
%
%% Plot the regridded variable
%figure(4)
%plotGrid(LLData4,'zonal');
%title('Data regridded from c24 to 4x5, Zonal')
%
%% Extract one layer of the plot
%figure(5)
%plotGrid(LLData4(:,:,1),'layer');
%title('Data regridded from c24 to 4x5, Level 1')
%
%% Compare lat-lon vs CS framework grid areas (compare sums!)
%gSpec = genGridSpec('gmao4x5','geos5'); 
%ll4_total_area = sum(gSpec.gridArea(:))
%pct_diff_cs_vs_4x5 = (ll4_total_area - cs_total_area) ...
%                     / ll4_total_area * 100




