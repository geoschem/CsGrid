function [ dxGrid,dyGrid ] = calcDXDY( lonEdge,latEdge )
%CALCDXDY Summary of this function goes here
%   Detailed explanation goes here
nP = size(lonEdge,1) - 1;
dxGrid = zeros(nP,nP,6);
dyGrid = zeros(nP,nP,6);
% Using FV3 rEarth for now (m)
rEarth = 6371e6;
for iTile = 1:6
    lonTile = lonEdge(:,:,iTile);
    latTile = latEdge(:,:,iTile);
    for iX = 1:nP
        for iY = 1:nP
            LL1 = [lonTile(iX,iY),latTile(iX,iY)];
            LL2 = [lonTile(iX+1,iY),latTile(iX+1,iY)];
            dxGrid(iX,iY,iTile) = greatCircleDistance(LL2,LL1,rEarth);
            LL2 = [lonTile(iX,iY+1),latTile(iX,iY+1)];
            dyGrid(iX,iY,iTile) = greatCircleDistance(LL2,LL1,rEarth);
        end
    end
end

end

