function gridArea = GEOSGridArea(latEdge,lonEdge)

if ischar(latEdge)
    % Try parsin it
    [lonEdge,latEdge] = parseGrid(latEdge);
end

% Return grid areas in m2
rEarth = 6.375e6; % in m
nLat = length(latEdge) - 1;
nLon = length(lonEdge) - 1;
%cellSAConst = 2*pi*rEarth*rEarth/nLon;
cellSAConst = rEarth * rEarth * pi * (lonEdge(end)-lonEdge(1)) / 180;
gridArea = ones(nLon,nLat);
for iLat = 1:nLat
    sinDiff = sind(latEdge(iLat+1))-sind(latEdge(iLat));
    cellSA = sinDiff.*cellSAConst; % in m2
    gridArea(:,iLat) = gridArea(:,iLat) .* cellSA;
end
end

