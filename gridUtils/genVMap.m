function [ hAxes, hPlot ] = genVMap( altEdge, latVals, varargin)
%GENVMAP Generate zonal axes with empty data
%   Inputs:
%       altEdge:    Vector of longitudes (degrees)
%       latVals:    Vector of latitudes (degrees)
%       Variable:
%           'position':     Axes position limits ([0.1 .1 .8 .8])
%           'fighandle':    Figure handle (gcf)
%           'subrange':     Altitude/latitude limits ([0 80],[-90 90])
%   Outputs:
%       hAxes:      Axes handle
%       hPlot:      Plot handle

% Parse input arguments
nArg = length(varargin);
iArg = 1;

% Set default values
axPos = [.1 .1 .8 .8];
altLims = [min(altEdge) max(altEdge)];
latLims = [min(latVals) max(latVals)];
hFig = gcf;
bigFont = 16;
smallFont = 12;
gridCells = true;

while iArg <= nArg
    currArg = varargin{iArg};
    switch lower(currArg)
        case 'position'
            axPos = varargin{iArg+1};
            nPlus = 2;
        case 'subrange'
            altLims = varargin{iArg+1};
            latLims = varargin{iArg+2};
            nPlus = 3;
        case 'fighandle'
            newFig = varargin{iArg+1};
            if ishandle(newFig)
                hFig = newFig;
            else
                error('genVAxes:badHandle','Specific figure handle does not exist.');
            end
            nPlus = 2;
        case 'smallfont'
            smallFont = varargin{iArg+1};
            nPlus = 2;
        case 'bigfont'
            bigFont = varargin{iArg+1};
            nPlus = 2;
        case 'gridcells'
            gridCells = varargin{iArg+1};
            nPlus = 2;
        otherwise
            warning('genVAxes:badArg','Argument ''%s'' not recognized. Ignoring.',...
                currArg);
            nPlus = 1;
    end
    iArg = iArg + nPlus;
end

% Set active figure
hAxes = axes('parent',hFig,'position',axPos);

% Generate raster reference
nLat = length(latVals);
nLev = length(altEdge);

dummyData = zeros(nLev,nLat);
dummyData(1,1) = 1;

if gridCells
    hPlot = pcolor(hAxes,latVals,altEdge,dummyData);
else
    [~,hPlot] = contourf(hAxes,latVals,altEdge,dummyData);
end
set(hPlot,'edgecolor','none');

set(hAxes,'layer','top','XGrid','on','YGrid','on',...
    'XLim',latLims,'YLim',altLims,'box','on','fontsize',smallFont);
xlabel(hAxes,'Latitude','interpreter','latex','fontsize',bigFont);
ylabel(hAxes,'Altitude (km)','interpreter','latex','fontsize',bigFont);

end