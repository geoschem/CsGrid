function [ gArea ] = calcGridArea( lonEdge, latEdge )
%CALCGRIDAREA Calculates grid areas (m2) for a rectangular grid

% Edge vectors are expected to be in degrees
rEarth = 6.375e6; % in m
nLat = length(latEdge) - 1;
nLon = length(lonEdge) - 1;
cellSAConst = 2*pi*rEarth*rEarth/nLon;
gArea = ones(nLon,nLat);
for iLat = 1:nLat
    sinDiff = sind(latEdge(iLat+1))-sind(latEdge(iLat));
    cellSA = sinDiff.*cellSAConst; % in m2
    gArea(:,iLat) = gArea(:,iLat) .* cellSA;
end

end

