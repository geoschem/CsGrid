function [ ] = setMapData( hPlot, gridData, varargin )
%SETMAPDATA Apply 2-D data to a horizontal or vertical grid
%   Inputs:
%       hPlot:      Handle for the plot
%       gridData:   MxN plot data. Either:
%                       [Lat x Lon]
%                       [Alt x Lat]

% Set default arguments
setAlt = false;

nArg = length(varargin);
iArg = 1;
reshapeLon = false;
while iArg <= nArg
    currArg = varargin{iArg};
    switch lower(currArg)
        case 'setz'
            % Update Z data as well
            setAlt = true;
            nPlus = 1;
        case 'lonshape'
            % Need to apply a longitudinal reshaping vector
            lonShapeVec = varargin{iArg+1};
            reshapeLon = length(lonShapeVec)>1;
            nPlus = 2;
        otherwise
            warning('setMapData:badArg','Argument ''%s'' not recognized.',...
                currArg);
            nPlus = 1;
    end
    iArg = iArg + nPlus;
end

% Apply longitudinal reshaping if required
if reshapeLon
    nLat = size(gridData,1);
    nLon = length(lonShapeVec);
    newData = zeros(nLat,nLon);
    for iLat = 1:nLat
        %miniVec = gridData(iLat,:);
        %newData(iLat,:) = miniVec(lonShapeVec);
        newData(iLat,:) = gridData(iLat,lonShapeVec);
    end
    gridData = newData;
end

% Pad grid data if necessary
nRowData = size(gridData,1);
nColData = size(gridData,2);

targField = 'CData';
if ischar(get(hPlot,targField))
    targField = 'ZData';
end
oldSize = size(get(hPlot,targField));
nRowPlot = oldSize(1);
nColPlot = oldSize(2);

sameSize = (nRowData == nRowPlot) && (nColData == nColPlot);
plusSize = (nRowData == (nRowPlot-1)) && (nColData == (nColPlot-1));

if plusSize
    % Determine a good mid-value
    dataRange = [min(gridData(:)),max(gridData(:))];
    if (dataRange(2) > 0) && (dataRange(1) < 0)
        edgeVal = 0;
    else
        edgeVal = mean(dataRange);
    end
    newData = edgeVal.*ones(oldSize);
    newData(1:nRowData,1:nColData) = gridData;
    gridData = newData;
    clear newData;
elseif ~sameSize
    % Data does not conform
    error('setMapData:wrongSize','Input data not sized correctly.');
end

% Set color data
set(hPlot,targField,gridData);

if setAlt
    % Set ZData as well
    set(hPlot,'ZData',gridData-max(gridData(:)));
end

end

