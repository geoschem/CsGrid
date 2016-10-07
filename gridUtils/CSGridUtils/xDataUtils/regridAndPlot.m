function [ dataOut ] = regridAndPlot( dataIn, xDataIn, convDir, ...
				      figTitle, varargin )
% SWAPCSFACES Transpose xData struct cubed sphere indexes for each face in
%             a vector of user-specified faces.
%
%    Input:  (1) dataIn: input data, corresponding to first 2 characters
%                        of convDir (e.g. cubed sphere if CS2LL)
%            (2) xDataIn: xData structure array for the regridding (see note)
%            (3) convDir: 'LL2CS' or 'CS2LL'
%            (4) figTitle: title for figure
%            (5) optional logical argument, turns off regridding if false
%
%    Output: (1) dataOut: conservatively regridded data, or the input
%                         data if regridding is turned off
%
%    NOTES: 
%        (1) xDataIn must have lat/lon info in first struct, xDataIn(1),
%            and cubed sphere info in second struct, xDataIn(2).
%        (2) If regridding is skipped by passing optional 5th arg false,
%            then assumes regridding has already occurred. Thus the input
%            data must have dimensions of the expected regridded data,
%            e.g. CS if passed 'LL2CS', or lat/lon if passed 'CS2LL'
%
% Lizzie Lundgren, 10/6/16

% allow user to turn off regridding with optional 5th argument (logical)
if nargin == 5
   regrid = varargin{1};
else
   regrid = true;  % regrid is true by default (on)
end

% Conservatively regrid. If regridding is skipped by passing false,
% assumes input data is already regridded, e.g. in CS if passed 'LL2CS'
if regrid
  [ dataOut, xDataOut ] = regridConservative( dataIn , xDataIn );
else
  dataOut = dataIn;
end

% Set variables for plotting face
numFaces = 6;
if strcmp(convDir,'LL2CS')
  NX = xDataIn(2).NX;
  xlim_max = NX;
  ylim_max = NX;
elseif strcmp(convDir,'CS2LL')
  xlim_max = xDataIn(1).NX;
  ylim_max = xDataIn(1).NY;
else
  error('Invalid conversion passed to regridAndPlot.m')
end
cbar_min = min(dataOut(:));
cbar_max = max(dataOut(:));

% Plot
figure;
for i=1:numFaces
    if strcmp(convDir,'LL2CS')
      data_to_plot = dataOut(:,i*NX-NX+1:i*NX)'; 
      subplot(2,3,i)
    else
      data_to_plot = dataOut';
    end    
    surf(data_to_plot,'EdgeColor','None');
    colorbar
    caxis([cbar_min,cbar_max])
    xlim([1,xlim_max])
    ylim([1,ylim_max])
    view(2); 
    title(sprintf('panel %d',i));
end   
suptitle(figTitle)
