function [ ] = plotLatLon( latEdge,lonEdge,Z,newPlot )
%PLOTLATLON Plot lat/lon data, given the edge vectors

if nargin > 3 && ~isempty(newPlot)
    assert(ischar(newPlot),'plotLatLon:badArg','newPlot argument must be a string or empty');
    clf;
    if strcmpi(newPlot,'flat')
        newPlot = 'robinson';
    end
    is3D = strcmpi(newPlot,'globe');
    axesm('globe','geoid',3670e3);
else
    currMap = gcm;
    mapAx = gca;
    is3D = strcmp(currMap.mapprojection,'globe');
    if ~is3D
        setm(mapAx,'MapProjection','globe');
    end
end

lonEdge = mod(lonEdge + 180,360) - 180;

% How many faces need to be split?
nSplit = sum(lonEdge < -180) * (length(latEdge) - 1);
nFacesBig = numel(Z) + nSplit;
latBig = zeros(4,nFacesBig);
lonBig = zeros(4,nFacesBig);
dataBig = zeros(nFacesBig,1);

hiFace = 0;
for iLon = 1:size(Z,1)
    iFace = 0;
    loLon = lonEdge(iLon);
    hiLon = lonEdge(iLon+1);
    lonSplit = loLon > hiLon;
    nFaces = size(Z,2) * (lonSplit + 1);
    latVec = zeros(4,nFaces);
    lonVec = zeros(4,nFaces);
    dataVec = zeros(nFaces,1);
    if lonSplit
        loLon = loLon - 360;
    end
    lonCorner = [loLon;loLon;hiLon;hiLon];
    for iLat = 1:size(Z,2)
        iFace = iFace + 1;
        latCorner = zeros(4,1);
        for iVert = 1:4
            xLat = iLat + (iVert == 1 || iVert == 4);
            latCorner(iVert) = latEdge(xLat);
        end
        lonVec(:,iFace) = lonCorner;
        latVec(:,iFace) = latCorner;
        cellVal = Z(iLon,iLat);
        dataVec(iFace) = cellVal;
        if lonSplit
            iFace = iFace + 1;
            lonVec(:,iFace) = lonCorner + 360;
            latVec(:,iFace) = latCorner;
            dataVec(iFace) = cellVal;
        end
    end
    loFace = hiFace + 1;
    hiFace = hiFace + iFace;
    latBig(:,loFace:hiFace) = latVec;
    lonBig(:,loFace:hiFace) = lonVec;
    dataBig(loFace:hiFace) = dataVec;
end

patchm(latBig,lonBig,'facevertexcdata',dataBig,...
        'EdgeColor','none','facecolor','flat');
    
maxZVal = 0;
if ~is3D
    mapAx = gca;
    setm(mapAx,'MapProjection','robinson');
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
    mlabel on;
    plabel on;
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
end

persistent coasts
if isempty(coasts)
    coasts = load('coast');
end
hold on;
plotm(coasts.lat,coasts.long,ones(size(coasts.lat)).*maxZVal,'-k');
hold off;
    
end

