function [ patchVec, mapAx ] = plotCSLayer( CSData, varargin )
%PLOTCSLAYER Plots CS data onto map axes (horizontal layer)

global CSGridDir

% Is the CS data separated into faces?
CSDims = size(CSData);
if length(CSDims) == 3
    assert(CSDims(2) == CSDims(1),'plotCSLayer:nonSquareSide','3D CS data must have square faces (total size NxNx6)');
    nPerSide = CSDims(1);
    reorderDefault = false;
elseif length(CSDims) == 2
    assert(CSDims(2) == 6*CSDims(1),'plotCSLayer:nonSquareSide','2D CS data must have square faces (total size Nx6N)');
    nPerSide = CSDims(1);
    reorderDefault = true;
else
    error('plotCSLayer:badDataDims','CS data must be 3D ([N x N x 6]) or 2D ([N x 6N]) in size');
end

% Generate the grid, if not supplied
inParse = inputParser;
plotOpts = {'mid','edge'};
iMid = 1;
iEdge = 2;
addParameter(inParse,'RebuildVtxData',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
addParameter(inParse,'QuickPlot',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
addParameter(inParse,'latEdge3D',[],@(x)validateattributes(x,{'numeric'},{'size',[nPerSide,nPerSide,6]}));
addParameter(inParse,'latEdge2D',[],@(x)validateattributes(x,{'numeric'},{'size',[nPerSide,nPerSide*6]}));
addParameter(inParse,'lonEdge3D',[],@(x)validateattributes(x,{'numeric'},{'size',[nPerSide,nPerSide,6]}));
addParameter(inParse,'lonEdge2D',[],@(x)validateattributes(x,{'numeric'},{'size',[nPerSide,nPerSide*6]}));
addParameter(inParse,'showEdge',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
addParameter(inParse,'parent',[],@ishandle);
addParameter(inParse,'offsetCube',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
addParameter(inParse,'projection','globe',@ischar);
addParameter(inParse,'style','edge',@(x)any(strcmpi(x,plotOpts)));
addParameter(inParse,'faceOrder',1:6,@(x)validateattributes(x,{'numeric'},{'integer','numel',6,'vector','positive','<',7}));
addParameter(inParse,'faceRotation',zeros(1,6),@(x)validateattributes(x,{'numeric'},{'integer','numel',6,'vector'}));
addParameter(inParse,'reorderData',reorderDefault,@(x)validateattributes(x,{'logical'},{'scalar'}));
addParameter(inParse,'coasts',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
addParameter(inParse,'labels',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
addParameter(inParse,'verbose',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
addParameter(inParse,'Renderer','opengl',@(x)validateattributes(x,{'char'},{}));
addParameter(inParse,'MapArg',{},@(x)validateattributes(x,{'cell'},{}));
parse(inParse,varargin{:});
offsetCube = inParse.Results.offsetCube;
verbose = inParse.Results.verbose;
mapArg = inParse.Results.MapArg;

plotStyleStr = inParse.Results.style;
plotStyle = find(strcmpi(plotStyleStr,plotOpts));
assert(~isempty(plotStyle),'plotCSLayer:badStyle','Unknown style ''%s'' requested',plotStyleStr);

flatPlot = false;
mapAx = inParse.Results.parent;
mapFound = false;
axFound = false;
if ~isempty(mapAx)
    % Is it a map axes?
    if ismap(mapAx)
        mapFound = true;
    elseif ishandle(mapAx)
        % Assume it's a standard set of axes
        axFound = true;
    else
        error('plotCSLayer:invalidParent','Parent must be map axes or standard axes');
    end
end

if ~mapFound
    if verbose
        fprintf('Generating map axes... ');
    end
    if ~axFound
        axArgs = {};
        clf;
    else
        % Inherit position from the given axes
        axArgs = {'position',get(mapAx,'position')};
    end
    set(gcf,'renderer',inParse.Results.Renderer);
    switch lower(inParse.Results.projection)
        case 'globe'
            % Project onto globe
            mapAx = axesm('globe','geoid',3670e3,axArgs{:});
        otherwise
        %case 'flat'
            % Roll the dice - assume it's a flat projection
            % If we're using midpoints, plot straight onto a flat
            % projection. If we're using edges, plot onto a globe and then
            % convert axes later
            targProj = lower(inParse.Results.projection);
            if strcmpi(targProj,'flat')
                targProj = 'Robinson';
            end
            if plotStyle == iEdge
                mapAx = axesm('globe','geoid',3670e3);
                if ~isempty(axArgs)
                    set(mapAx,axArgs{:});
                end
            else
                % Generate a horizontal projection
                if axFound
                    axes(mapAx);
                end
                mapAx = worldmap('world',axArgs{:});
            end
            flatPlot = true;
        %otherwise
            %error('plotCSLayer:unknownProjection','Projection not recognized');
    end
    if verbose
        fprintf('done.\n');
    end
end

lat3D = ~isempty(inParse.Results.latEdge3D);
lat2D = ~isempty(inParse.Results.latEdge2D);
lon3D = ~isempty(inParse.Results.lonEdge3D);
lon2D = ~isempty(inParse.Results.lonEdge2D);

edgeOptSum = lat3D+lat2D+lon3D+lon2D;
switch plotStyle
    case iMid
        if edgeOptSum ~= 0
            warning('plotCSLayer:automidOnly','Midpoint plotting must be performed with auto-grid generation - ignoring user grid');
            lat3D = false;
            lon3D = false;
            lat2D = false;
            lon2D = false;
            edgeOptSum = 0;
        end
        [~,~,lonMid,latMid] = calcCSGrid(nPerSide,'offsetcube',offsetCube);
end

% Do we need to convert 2D to 3D data?
%{
faceOrder = inParse.Results.faceOrder;
faceRot = inParse.Results.faceRotation;
if GCHPData
    % Override
    faceOrder = [1,2,6,4,5,3];
    faceRot = [0,0,90,180,180,90];
end
assert(all(sort(faceOrder)==(1:6)),'plotCSLayer:badFaceOrder','Face order must contain all integers from 1 to 6');
if grid2Dto3D
    CS2D = CSData;
    CSData = zeros(nPerSide,nPerSide,6);
    for iFace = 1:6
        % Change READ order
        readFace = faceOrder(iFace);
        faceLo = (nPerSide * (readFace-1)) + 1;
        faceHi = faceLo + nPerSide - 1;
        temp2D = CS2D(:,faceLo:faceHi);
        currRot = faceRot(iFace);
        n90 = round(currRot/90);
        if n90 ~= 0
            temp2D = rot90(temp2D,n90);
        end
        CSData(:,:,iFace) = temp2D;
    end
    clear CS2D;
end
%}
if inParse.Results.reorderData
    CSData = reorderCS(CSData,inParse.Results.reorderData);
end

assert(any(edgeOptSum == [0,2]),'plotCSLayer:badEdgeSpec','Invalid edge specification for grid');
if edgeOptSum == 0
    % Generate grid
    [lonEdge,latEdge] = calcCSGrid(nPerSide,'offsetcube',offsetCube);
elseif edgeOptSum == 2
    if (lon3D && lat3D)
        lonEdge = inParse.Results.lonEdge3D;
        latEdge = inParse.Results.latEdge3D;
    elseif (lon2D && lat2D)
        lonEdge2D = inParse.Results.lonEdge2D;
        latEdge2D = inParse.Results.latEdge2D;
        nP1 = nPerSide + 1;
        lonEdge = zeros(nP1,nP1,6);
        latEdge = zeros(nP1,nP1,6);
        faceHi = 0;
        for iFace = 1:6
            faceLo = faceHi + 1;
            faceHi = faceHi + nP1;
            lonEdge(:,:,iFace) = lonEdge2D(:,faceLo:faceHi);
            latEdge(:,:,iFace) = latEdge2D(:,faceLo:faceHi);
        end
        clear lonEdge2D latEdge2D;
    else
        error('plotCSLayer:mismatchEdgeSpec','Please specify either entirely 2D edges, entirely 3D edges, or none');
    end
end

assert(all(size(lonEdge)==size(latEdge)),'plotCSLayer:mistmatchLonLat','Longitude and latitude edges mismatched');
assert(all(size(lonEdge)==([1,1,0]+size(CSData))),'plotCSLayer:dataEdgeMismatch','CS data does not conform to specified grid');

% Ensure that there is something to plot
cLim = [min(CSData(:)),max(CSData(:))];
if diff(cLim)==0
    miniDiff = [-1.0e-15,1.0e-15];
    if abs(cLim(1)) <= miniDiff(2)
        cLim = miniDiff;
    else
        cLim = cLim(1).*(1.0+miniDiff);
    end
end

set(mapAx,'clim',cLim);

% Skip patch validity check
quickPlot = inParse.Results.QuickPlot;

% For later updates
patchVec = zeros(6,1);
if verbose
    fprintf('Plotting faces: [');
end
if plotStyle == iMid
    % Do each face separately - anything else causes.. issues
    for iFace = 1:6
        if verbose
            if iFace < 6
                fprintf('%i,',iFace);
            else
                fprintf('%i]... ',iFace);
            end
        end
        lonGrid = lonMid(:,:,iFace);
        latGrid = latMid(:,:,iFace);
        dataGrid = CSData(:,:,iFace);
        patchVec(iFace) = geoshow(latGrid,lonGrid,dataGrid,...
            'displaytype','texturemap');
    end
elseif plotStyle == iEdge
    % Do we want to show edges?
    if inParse.Results.showEdge
        edgeColor = [1 1 1].*0.5;
    else
        edgeColor = 'none';
    end
    % Draw each one as a patch using edge data
    % Do we have pre-calculated vertex data?
%    vtxDataFile = fullfile('GridData','VertexData',...
%        sprintf('VtxData_C%i.mat',nPerSide));
    vtxDataFile = fullfile(CSGridDir,'GridData','VertexData',...
        sprintf('VtxData_C%i.mat',nPerSide));
    % Allow for cube offset
    lonMatName = sprintf('lonMat_Offset%d',offsetCube);
    latMatName = sprintf('latMat_Offset%d',offsetCube);
    storedVertexData = fastExist(vtxDataFile) && ~inParse.Results.RebuildVtxData;
    % Check that the wanted files also exist
    if storedVertexData
        dataList = who('-file',vtxDataFile);
        storedVertexData = any(strcmpi(dataList,lonMatName)) && ...
            any(strcmpi(dataList,latMatName));
    end
    if storedVertexData
        % Read the vertex data
        fileData = load(vtxDataFile,lonMatName,latMatName);
        lonMat = fileData.(lonMatName);
        latMat = fileData.(latMatName);
        clear fileData;
    else
        % Save the data for storage
        lonMat = nan(4,nPerSide*nPerSide,6);
        latMat = nan(4,nPerSide*nPerSide,6);
    end
    for iFace = 1:6
        if verbose
            if iFace < 6
                fprintf('%i,',iFace);
            else
                fprintf('%i]... ',iFace);
            end
        end
        faceData = CSData(:,:,iFace);
        faceLon = lonEdge(:,:,iFace);
        faceLat = latEdge(:,:,iFace);
        dataVec = nan(nPerSide*nPerSide,1);
        if ~storedVertexData
            lonVec = nan(4,nPerSide*nPerSide);
            latVec = nan(4,nPerSide*nPerSide);
            iIndex = 0;
            for iLon = 1:nPerSide
                for iLat = 1:nPerSide
                    latCorner = zeros(4,1);
                    lonCorner = zeros(4,1);
                    for iVert = 1:4
                        xLon = iLon + (iVert > 2);
                        xLat = iLat + (iVert == 1 || iVert == 4);
                        lonCorner(iVert) = faceLon(xLon,xLat);
                        latCorner(iVert) = faceLat(xLon,xLat);
                    end
                    if flatPlot && ~quickPlot
                        [lonCorner,latCorner] = validLLPatch(lonCorner,latCorner);
                    end
                    iIndex = iIndex + 1;
                    lonVec(:,iIndex) = lonCorner;
                    latVec(:,iIndex) = latCorner;
                    dataVec(iIndex) = faceData(iLon,iLat);
                end
            end
            lonMat(:,:,iFace) = lonVec;
            latMat(:,:,iFace) = latVec;
        else
            lonVec = lonMat(:,:,iFace);
            latVec = latMat(:,:,iFace);
            iIndex = 0;
            for iLon = 1:nPerSide
                for iLat = 1:nPerSide
                    iIndex = iIndex + 1;
                    dataVec(iIndex) = faceData(iLon,iLat);
                end
            end
        end
        if flatPlot
            patchVec(iFace) = patchm(latVec,lonVec,0,'facevertexcdata',dataVec,...
                    'facecolor','flat','EdgeColor',edgeColor);
        else
            patchVec(iFace) = patchm(latVec,lonVec,'facevertexcdata',dataVec,...
                    'facecolor','flat','EdgeColor',edgeColor);
        end
    end
    if ~storedVertexData && ~quickPlot
        fileData.(latMatName) = latMat;
        fileData.(lonMatName) = lonMat;
        if fastExist(vtxDataFile)
            % Must not contain the data we need
            save(vtxDataFile,'-struct','fileData','-append');
        else
            save(vtxDataFile,'-struct','fileData');
        end
    end
end
if verbose
    fprintf('done.\n');
end

% Flatten the plot, if necessary
maxZVal = 1;
if flatPlot
    if verbose
        fprintf('Flattening map projection... ');
    end
    setm(mapAx,'MapProjection',targProj,mapArg{:});
%     for iArg = 1:length(mapArg)
%         disp(mapArg{iArg});
%     end
%     pause;
    mapChildren = get(gca,'Children');
    maxZVal = 0;
    for iChild = 1:length(mapChildren)
        if any(strcmpi(mapChildren(iChild).Type,{'Patch','Surface'}))
            maxZVal = max([max(mapChildren(iChild).ZData(:)),maxZVal]);
        end
    end
    maxZVal = maxZVal + 1;
    if verbose
        fprintf('done.\n');
    end
end
if flatPlot
    framem on;
    gridm on;
    % Do we want labels?
    if inParse.Results.labels
%        setm(gca,'parallellabel','off');
%        setm(gca,'meridianlabel','off');
        mlabel on;
        plabel on;
    end
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
    set(gca,'visible','off');
end

persistent coasts
if inParse.Results.coasts
    if verbose
        fprintf('Overlaying coasts... ');
    end
    if isempty(coasts)
        coasts = load('coast');
    end
    hold on;
    plotm(coasts.lat,coasts.long,ones(size(coasts.lat)).*maxZVal,'-k');
    hold off;
    if verbose
        fprintf('done.\n');
    end
end

end
