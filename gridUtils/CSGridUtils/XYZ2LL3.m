function [lonGrid,latGrid,xyzArr] = XYZ2LL3(xyzGrid)
%XYZ2LL3 Convert 2-D grid of cartesian to lat-lon co-ords
% This is just a wrapper function for XYZ2LL2, which handles a single point
nX = size(xyzGrid,2);
nY = size(xyzGrid,3);
lonGrid = zeros(nX,nY);
latGrid = zeros(nX,nY);
xyzArr = xyzGrid;
% Go column by column
for iY = 1:nY
    [lonGrid(:,iY),latGrid(:,iY),xyzArr(:,:,iY)] = XYZ2LL2(xyzGrid(:,:,iY));
end
end