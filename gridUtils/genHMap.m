function [ hAxes, hPlot,latEdge, lonEdge, lonShapeVec ] = genHMap( latEdge, lonEdge, varargin)
%GENHMAP Generate horizontal axes with empty data
%   Inputs:
%       latEdge:    Vector of latitudes (degrees)
%       lonEdge:    Vector of longitudes (degrees)
%       Variable:
%           'position':     Axes position limits ([0.1 .1 .8 .8])
%           'fighandle':    Figure handle (gcf)
%           'origin':       Map origin ([0 0])
%           'subrange':     Latitude/longitude limits ([-90 90],[-180 180])
%           'type':         Display type ('texturemap')
%   Outputs:
%       hAxes:      Axes handle
%       hPlot:      Plot handle
%       lonShapeVec:Vector mapping original longitudes to new ones, if req.

% Parse input arguments
nArg = length(varargin);
iArg = 1;

% Set default values
%posArg = {};
projID = 'wagner4';
latLims = [min(latEdge) max(latEdge)];
lonLims = [min(lonEdge) max(lonEdge)];
mapOrigin = [0 0];
smallFont = 12;
MLineLoc = 30;
PLineLoc = 15;
noLabel = true;
showCoast = true;
coastColor = [0,0,0];

while iArg <= nArg
    currArg = varargin{iArg};
    switch lower(currArg)
        case 'coastcolor'
            coastColor = varargin{iArg+1};
            nPlus = 2;
        case 'nocoast'
            showCoast = false;
            nPlus = 1;
        case 'mlineloc'
            MLineLoc = varargin{iArg+1};
            nPlus = 2;
        case 'plineloc'
            PLineLoc = varargin{iArg+1};
            nPlus = 2;
        case 'position'
            axPos = varargin{iArg+1};
            posArg = {'position',axPos};
            nPlus = 2;
        case 'projection'
            projID = varargin{iArg+1};
            nPlus = 2;
        case 'origin'
            mapOrigin = varargin{iArg+1};
            nPlus = 2;
        case 'subrange'
            lonLims = varargin{iArg+1};
            latLims = varargin{iArg+2};
            nPlus = 3;
        case 'fighandle'
            newFig = varargin{iArg+1};
            if ishandle(newFig)
                hFig = newFig;
            else
                error('genHAxes:badHandle','Specific figure handle does not exist.');
            end
            nPlus = 2;
        case 'axhandle'
            % Generate the map axes in the given axes
            testAx = varargin{iArg+1};
            if ishandle(testAx)
                hAxes = testAx;
            else
                nSkip = nSkip + 1;
            end
            nPlus = 2;
        case 'smallfont'
            smallFont = varargin{iArg+1};
            nPlus = 2;
        case 'bigfont'
            % Ignore
            nPlus = 2;
        case 'nolabel'
            noLabel = true;
            nPlus = 1;
        otherwise
            warning('genHAxes:badArg','Argument ''%s'' not recognized. Ignoring.',...
                currArg);
            nPlus = 1;
    end
    iArg = iArg + nPlus;
end

% Set active figure
if ~exist('hFig','var')
    if exist('hAxes','var')
        hFig = get(hAxes,'parent');
    else
        hFig = gcf;
    end
end
figure(hFig);

% Handle special case - full grid with looping data
lonShapeVec = 0;
lonMax = 180;
lonMin = -180;
eastExt = (lonLims(2)>lonMax);
westExt = (lonLims(1)<lonMin);
if westExt || eastExt
    % Produce new longitude vector
    if eastExt && westExt
        error('genHMap:badLims','Longitude values exceed 360 degree rotation.');
    elseif eastExt
        % Extra cell needed at western limit for each edge exceeding 180
        nExtra = sum(lonEdge>lonMax);
        nEdge = length(lonEdge);
        newLonVec = zeros(nEdge + 1,1);

        % First entry: -180
        newLonVec(1) = lonMin;%-180;
        % Second entry: First value to exceed 180, minus 360
        newLonVec(2:(nExtra+1)) = lonEdge(lonEdge>lonMax) - 360;
        newLonVec((nExtra+1):end-1) = lonEdge(lonEdge<=lonMax);
        newLonVec(end) = lonMax;
        
        %newLonVec(1) = newLonVec(2) - diff(newLonVec(2:3));
        %newLonVec(end) = newLonVec(end-1) + diff(newLonVec(end-2:end-1));
        
        % Vector will describe how to map old values to new
        nCell = nEdge - 1;
        % Mapping nCell values to nCell + 1, duplicating onces
        lonShapeVec = zeros(nCell + 1,1);
        % Duplicated value
        lonShapeVec(1:nExtra) = (nCell+1-nExtra):nCell;
        lonShapeVec(nExtra+1:end) = 1:(nCell+1-nExtra);
        
        lonEdge = newLonVec;
    elseif westExt
        % As for East, but reversed
        nExtra = sum(lonEdge<lonMin);
        nEdge = length(lonEdge);
        newLonVec = zeros(nEdge + 1,1);

        % Last entry: 180
        newLonVec(end) = lonMax;
        % Previous entries: those that are below 180, plus 360
        newLonVec(end-(1:nExtra)) = lonEdge(lonEdge<lonMin) + 360;
        newLonVec(2:(end-nExtra)) = lonEdge(lonEdge>lonMin);
        newLonVec(1) = lonMin;
        
        %newLonVec(1) = newLonVec(2) - diff(newLonVec(2:3));
        %newLonVec(end) = newLonVec(end-1) + diff(newLonVec(end-2:end-1));
        
        % Vector will describe how to map old values to new
        nCell = nEdge - 1;
        % Mapping nCell values to nCell + 1, duplicating onces
        lonShapeVec = zeros(nCell + 1,1);
        % Duplicated value
        lonShapeVec((end+1-nExtra):end) = 1:nExtra;
        lonShapeVec(1:(end-nExtra)) = nExtra:nCell;
        
        lonEdge = newLonVec;
    end
    lonLims = [min(lonEdge) max(lonEdge)];
end

% Generate a horizontal projection
hAxes = worldmap('world');
setm(hAxes,'MapProjection',projID);
mapChildren = get(gca,'Children');
maxZVal = 0;
for iChild = 1:length(mapChildren)
    if any(strcmpi(mapChildren(iChild).Type,{'Patch','Surface'}))
        maxZVal = max([max(mapChildren(iChild).ZData(:)),maxZVal]);
    end
end
maxZVal = maxZVal + 1;

framem on;
gridm on;

setm(gca,'MapLatLimit',latLims,'MapLonLim',lonLims,...
    'origin',mapOrigin,'GAltitude',0,'fontsize',smallFont,...
    'MLineLocation',MLineLoc,'PLineLocation',PLineLoc);

mapChildren = get(gca,'Children');
for iChild = 1:length(mapChildren)
    isFrame = strcmpi(mapChildren(iChild).Tag,'Frame');
    isLine = any(strcmpi(mapChildren(iChild).Tag,{'Parallel','Meridian'}));
    if isFrame || isLine
        oldData = mapChildren(iChild).ZData;
        set(mapChildren(iChild),'ZData',ones(size(oldData)) .* maxZVal);
        if isLine
            set(mapChildren(iChild),'color',[0 0 0]);
        end
    end
end

% Do we want labels?
if noLabel
    mlabel off;
    plabel off;
else
    mlabel on;
    plabel on;
end

%{
hAxes = axesm(projID,'MeridianLabel',strLabel,...
    'ParallelLabel',strLabel,'MapLatLimit',latLims,'MapLonLim',lonLims,...
    'Grid','on','frame','on','origin',mapOrigin,'GAltitude',0,...
    'fontsize',smallFont,'MLineLocation',MLineLoc,'PLineLocation',PLineLoc);
    %'fontsize',smallFont,'MLineLocation',90,'PLineLocation',45);
%}

% Generate raster reference
nLat = length(latEdge)-1;
nLon = length(lonEdge)-1;

dummyData = zeros(nLat,nLon);
dummyData(1,1) = 1;

% Force the plot data to be slightly below the coast data
hPlot = pcolorm(latEdge(1:nLat+1),lonEdge(1:nLon+1),dummyData,...
    'zdata',-0.1.*ones(nLat+1,nLon+1),'parent',hAxes);
uistack(hPlot,'bottom');
if showCoast
    landShapes = shaperead('landareas','UseGeoCoords',true);
    geoshow(hAxes,[landShapes.Lat], [landShapes.Lon], 'Color', coastColor);%,'linewidth',1);
end

set(gca,'visible','off');

% Setting position later seems to work better
%set(hAxes,'Visible','off','box','off',posArg{:});

end