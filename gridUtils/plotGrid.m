function [ hFig, hAxVec, hPlotVec ] = plotGrid( gridData, plotType, varargin )
%PLOTGRID High-level function to plot data
%   Inputs:
%       gridData:   Regular 3-D gridded data
%       plotType:   Either:
%                       'zonal':    Zonal mean
%                       'layer':    Horizontal model layer
%   At least one set of co-ordinates must also be supplied (altitude and
%   latitude for zonal plots, longitude and latitude for layer plots)
%   through the optional arguments. Alternatively, a standard resolution
%   can be applied through use of the 'res' option.
%   Optional arguments:
%       latEdge:    Vector of latitude edge values (degrees)
%       lonEdge:    Vector of longitude edge values (degrees)
%       altEdge:    Vector of altitude edge values (degrees)
%       

% Set up defaults
nArg = length(varargin);
iArg = 1;
iLayer = 1;
%plotUnit = '';
limSet = false;
smallFont = 12;
bigFont = 16;
subRangeSpec = false;
originSpec = false;
% 2016-04-01: Default to Robinson
customProj = true;
mapProj = 'Robinson';
limLev = false;
gridCells = true;
figName = [];
usePOpts = false;
nSkip = 0;
PLineLoc = 15;
MLineLoc = 30;
noLabel = true;
showContour = false;
noCoast = false;
gSpec = [];
hrzSpec = '';
vrtSpec = '';
coastColor = [0,0,0];

% SDE 2016-03-08
%   OpenGL is sometimes slightly slower to plot, but it places much less of
%   a burden on the system long-term - hugely preferable
newRenderer = 'opengl';

while iArg <= nArg
    currArg = varargin{iArg};
    if ~ischar(currArg)
        nPlus = 1;
        nSkip = nSkip + 1;
    else
        switch lower(currArg)
            case 'renderer'
                newRenderer = varargin{iArg+1};
                nPlus = 2;
            case 'coastcolor'
                coastColor = varargin{iArg+1};
                nPlus = 2;
            case 'nocoast'
                noCoast = true;
                nPlus = 1;
            case 'showcontour'
                % Overlay contour map
                showContour = true;
                contourArg = varargin{iArg+1};
                nPlus = 2;
            case 'transpose'
                gridData = transpose(gridData);
                nPlus = 1;
            case 'showlabels'
                noLabel = false;
                nPlus = 1;
            case 'mlineloc'
                MLineLoc = varargin{iArg+1};
                nPlus = 2;
            case 'plineloc'
                PLineLoc = varargin{iArg+1};
                nPlus = 2;
            case 'plotopts'
                plotOpts = varargin{iArg+1};
                bigFont = plotOpts.bigFont;
                smallFont = plotOpts.smallFont;
                usePOpts = true;
            case {'newfig','forcenewfig'}
                hFig = figure('color',[1 1 1],'windowstyle','docked');
                nPlus = 1;
            case 'gridcells'
                gridCells = varargin{iArg + 1};
                nPlus = 2;
            case {'targlayer','ilayer','layer'}
                iLayer = varargin{iArg+1};
                nPlus = 2;
            case 'squarecell'
                cellSize = varargin{iArg+1};
                %{
                nLat = round(180/cellSize);
                nLon = 2*nLat;
                latEdge = linspace(-90,90,nLat+1);
                lonEdge = linspace(-180,180,nLon+1);
                %}
                hrzSpec = sprintf('SQR_%i',cellSize);
                nPlus = 2;
            case 'cellsize'
                % [ lonSize latSize ]
                cellSize = varargin{iArg+1};
                lonCellSize = cellSize(1);
                latCellSize = cellSize(2);
                %{
                nLon = round(360/lonCellSize);
                nLat = round(180/latCellSize);
                lonEdge = linspace(-180,180,nLon+1);
                latEdge = linspace(-90,90,nLat+1);
                %}
                hrzSpec = sprintf('RECT_%ix%i',latCellSize,lonCellSize);
                nPlus = 2;
            case 'gridspec'
                gSpec = varargin{iArg+1};
                if isa(gSpec,'modelData')
                    gSpec = gSpec.gSpec;
                elseif ~isa(gSpec,'gridSpec')
                    error('plotGrid:badGridSpec','Argument following GridSpec argument must be either a modelData or gridSpec object');
                end
                %{
                latEdge = gSpec.latEdge;
                lonEdge = gSpec.lonEdge;
                altEdge = gSpec.zEdge;
                %}
                nPlus = 2;
                %{
            case {'latedge','latvals'}
                latEdge = varargin{iArg+1};
                nPlus = 2;
            case {'altedge','altvals'}
                altEdge = varargin{iArg+1};
                nPlus = 2;
            case {'lonedge','lonvals'}
                lonEdge = varargin{iArg+1};
                nPlus = 2;
                %}
            case 'maxlev'
                limLev = true;
                maxLev = varargin{iArg+1};
                nPlus = 2;
            case {'subrange'}
                subLon = varargin{iArg+1};
                subLat = varargin{iArg+2};
                subRangeSpec = true;
                nPlus = 3;
            case {'origin'}
                originSpec = true;
                plotOrigin = varargin{iArg+1};
                nPlus = 2;
            case {'res'}
                %[lonEdge,latEdge,altEdge,isNested] = parseGrid(varargin{iArg+1},numel(gridData));
                hrzSpec = varargin{iArg+1};
                vrtSpec = varargin{iArg+2};
                %{
                gSpec = genGridSpec(hrzSpec,vrtSpec);
                latEdge = gSpec.latEdge;
                lonEdge = gSpec.lonEdge;
                altEdge = gSpec.zEdge;
                %}
                nPlus = 3;
            case 'bigfont'
                bigFont = varargin{iArg+1};
                nPlus = 2;
            case 'smallfont'
                smallFont = varargin{iArg+1};
                nPlus = 2;
                %{
            case {'unit','plotunit'}
                plotUnit = varargin{iArg+1};
                nPlus = 2;
            case 'nbin'
                nBin = varargin{iArg + 1};
                nPlus = 2;
                %}
            case 'clim'
                cLims = varargin{iArg+1};
                limSet = true;
                nPlus = 2;
            case 'fighandle'
                newFig = varargin{iArg+1};
                if ishandle(newFig)
                    hFig = newFig;
                else
                    error('plotGrid:badHandle','Specific figure handle does not exist.');
                end
                nPlus = 2;
            case 'axhandle'
                newAx = varargin{iArg+1};
                if ishandle(newAx)
                    hAx = newAx;
                    hFig = get(hAx,'parent');
                else
                    error('plotGrid:badHandle','Specific axis handle does not exist.');
                end
                nPlus = 2;
            case 'projection'
                customProj = true;
                mapProj = varargin{iArg+1};
                nPlus = 2;
            case 'name'
                figName = varargin{iArg+1};
                nPlus = 2;
            otherwise
                
                warning('plotGrid:badArg','Argument ''%s'' not recognized. Ignoring.',...
                    currArg);
                
                nPlus = 1;
                nSkip = nSkip + 1;
        end
    end
    iArg = iArg + nPlus;
end

if nSkip > 0
    warning('plotGrid:badArg','Skipped %i bad argument(s)',nSkip);
end

nDim = ndims(gridData);

% Get plot type
switch plotType
    case 'layer'
        isLayer = true;
        gridCells = true;
    case 'zonal'
        isLayer = false;
    otherwise
        error('plotGrid:badPlotType','Plot type ''%s'' not recognized.',...
            plotType);
end

assert(~(showContour&&isLayer),'plotGrid:contourLayer','Cannot show contour maps on layer plots');

% Do we have a grid specification?
if isempty(gSpec)
    if isLayer
        vrtSpec = [];
        if isempty(hrzSpec)
            hrzSpec = sprintf('auto_%ix%i',size(gridData,2),size(gridData,1));
        end
    else
        if isempty(hrzSpec)
            assert(nDim == 3,'plotGrid:cannotAutoCalculate',...
                'Cannot auto-determine grid specification from pre-averaged zonal data');
            hrzSpec = sprintf('auto_%ix%i',size(gridData,1),size(gridData,2));
        end
        if isempty(vrtSpec)
            vrtSpec = sprintf('auto_%i',size(gridData,nDim));
        end
    end
    [gSpec,doTranspose] = genGridSpec(hrzSpec,vrtSpec);
    if doTranspose
        gridData = gridData';
    end
end

% Retrieve grid properties from grid specification
latEdge = gSpec.latEdge;
lonEdge = gSpec.lonEdge;
altEdge = gSpec.zEdge;
isNested = gSpec.isNested;

% Nested grid?
if isNested && ~subRangeSpec
    % Make sure the output only shows the nested subrange
    subRangeSpec = true;
    subLon = lonEdge([1 end]);
	subLat = latEdge([1 end]);
end

% Only worry about altitude limits if not a layer plot
limLev = limLev && ~isLayer;

if nDim == 3
    if isLayer
        % Isolate specific layer
        gridData = transpose(gridData(:,:,iLayer));
	else
        % Take zonal mean
        gridData = transpose(squeeze(mean(gridData,1)));
    end
elseif nDim ~= 2
    error('plotGrid:badSize','Grid not correctly sized.');
end

% Cut off unwanted layers
if limLev
    gridData = gridData(1:maxLev,:);
    if gridCells
        altEdge = altEdge(1:maxLev+1);
    else
        altEdge = altEdge(1:maxLev);
    end
end

% Determine (or apply) plot limits
if limSet
    % Force grid into limits
    gridData(gridData<max(cLims)) = max(cLims);
    gridData(gridData>min(cLims)) = min(cLims);
else
    cLims = [min(gridData(:)) max(gridData(:))];
    if issparse(cLims)
        cLims = full(cLims);
    end
    
    % Safegard against NaN values
    gridNaN = isnan(gridData);
    gridInf = isinf(gridData);
    gridBad = gridNaN | gridInf;
    assert(~all(gridBad(:)),'plotGrid:noData','No plottable data');
    if any(gridBad(:))
        tempData = gridData(~gridBad);
        cLims(2) = max(tempData(:));
        cLims(1) = min(tempData(:));
    end
end

% Check grid dimensions
inSize = size(gridData);

if isLayer
    hEdge = lonEdge;
    vEdge = latEdge;
else
    hEdge = latEdge;
    vEdge = altEdge;
end

% Code gets upset if we use single precision
hEdge = cast(hEdge,'double');
vEdge = cast(vEdge,'double');

% Compensate for some tracers only corresponding to subranges
if (~isLayer) && length(vEdge) > inSize(1) + 1
    warning('plotGrid:badDim',...
        ['Fewer vertical levels than altitude edges. ',...
        'Assuming subrange from surface.']);
    vEdge = vEdge(1:(inSize(1)+1));
end

altMatch = inSize(1) + (isLayer||gridCells);
if (length(vEdge) ~= altMatch) || (length(hEdge) ~= (inSize(2)+1))
    % Does transpose work?
    if (length(hEdge) == altMatch) && (length(vEdge) == (inSize(2)+1))
        warning('plotGrid:dimTranspose','Using transposed grid');
        gridData = gridData';
    else
        error('plotGrid:badDim','Dimensions are not correctly set up.');
    end
end

% Sanitize the data - "infs" will cause a crash
infData = isinf(gridData);
if any(infData)
    warning('plotGrid:infData','Infinities found in data. Setting to NaN');
    gridData(infData) = nan;
end

% Set up figure
if ~exist('hFig','var')
    hFig = gcf;
    %set(hFig,'color',[1 1 1],'windowstyle','docked');
    set(hFig,'color',[1 1 1]);
end

if usePOpts
    plotOpts.addFig(hFig);
end

% Clear the current plot figure
extraArg = {};
if ~exist('hAx','var')
    clf(hFig);
else
    axPos = get(hAx,'position');
    extraArg = [extraArg {'position',axPos}];
    if isLayer
        extraArg = [extraArg {'axhandle',hAx}];
    else
        delete(hAx);
    end
end

% Set renderer just in case
set(hFig,'renderer',newRenderer);

% Plot grid data first

% Set up axes
%hAxVec = 0;
hPlotVec = 0;
if isLayer
    if customProj
        extraArg = [extraArg {'projection',mapProj}];
    end
    if subRangeSpec
        % Sub-range specified
        extraArg = [extraArg {'subrange',subLon,subLat}];
    end
    if originSpec
        extraArg = [extraArg {'origin',plotOrigin}];
    end
    if noLabel
        extraArg = [extraArg {'nolabel'}];
    end
    if noCoast
        extraArg = [extraArg {'nocoast'}];
    end
    [hAxVec,hPlotVec(1),~,~,lonShapeVec] = genHMap(vEdge,hEdge,...
        'fighandle',hFig,'smallfont',smallFont,'coastColor',coastColor,...
        'bigfont',bigFont,'PLineLoc',PLineLoc,'MLineLoc',MLineLoc,extraArg{:});
    mapArg = {'lonshape',lonShapeVec};
else
    if gridCells
        hVals = hEdge;
    else
        % Send hCenter, not hEdge
        hVals = hEdge(1:end-1) + (0.5.*diff(hEdge));
    end
    [hAxVec(1),hPlotVec(1)] = genVMap(vEdge,hVals,...
        'fighandle',hFig,'smallfont',smallFont,'bigfont',bigFont,...
        'gridCells',gridCells,extraArg{:});
    mapArg = {};
end

% Apply data to established axes
setMapData(hPlotVec(1),gridData,mapArg{:});

% Safeguard against uniform values
if (cLims(1) == cLims(2))
    cLims(1) = cLims(2) - 0.5;
    cLims(2) = cLims(2) + 0.5;
end

assert(cLims(2)>cLims(1) && ~any(isnan(cLims)|isinf(cLims)),'plotGrid:badCLim','Color limits must rise monotonically!');

set(hAxVec(1),'clim',cLims);

% Set colormap (jet, for now)
nCells = 401;
if verLessThan('matlab','8.3.0')
    cMap = jet(nCells);
else
    cMap = parula(nCells);
end
colormap(hAxVec(1),cMap);

if ~verLessThan('matlab','8.4.0')
    set(hAxVec(1),'ticklabelinterpreter','latex');
end

% Show a contour map if necessary
if showContour
    if ~iscell(contourArg)
        contourArg = {contourArg};
    end
    vMid = ((vEdge(1:(end-1))+vEdge(2:end))./2.0);
    hMid = ((hVals(1:(end-1))+hVals(2:end))./2.0);
    hold on;
    [conData,conHandle] = contour(hAxVec(1),hMid,vMid,gridData,contourArg{:},'-k');
    clabel(conData,conHandle);
    hold off;
end

% Set figure name, if requested
if ~isempty(figName)
    set(hFig,'name',figName);
end

end
