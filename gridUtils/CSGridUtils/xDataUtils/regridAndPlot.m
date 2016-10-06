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
%
%    Output: (1) dataOut: conservatively regridded data
%
%    NOTE: xDataIn must have lat/lon info in first struct, xDataIn(1),
%          and cubed sphere info in second struct, xDataIn(2).
%
% Lizzie Lundgren, 10/6/16

%% allow user to turn off regridding with optional 5th argument (logical)
% NOTE: this is not implemented yet
%if nargin == 5
%   regrid = varargin{1};
%else
%   regrid = true;  % regrid is true by default (on)
%end

% Conservatively regrid
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
    caxis([0,numFaces])
    xlim([1,xlim_max])
    ylim([1,ylim_max])
    view(2); 
    title(sprintf('panel %d',i));
end   
suptitle(figTitle)
