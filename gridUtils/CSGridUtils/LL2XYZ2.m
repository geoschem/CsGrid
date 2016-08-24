function [xyz] = LL2XYZ2(lonPt,latPt)
%LL2XYZ2 Converts a lon/lat pair (in radians) to cartesian co-ordinates
% Vector should point to the surface of the unit sphere
nPts = length(lonPt);
xyz = zeros(3,nPts);
for iPt = 1:nPts
    xyz(:,iPt) = LL2XYZ([lonPt(iPt),latPt(iPt)]);
end
end

function [xyz] = LL2XYZ(lonLatPair)
xPt = cos(lonLatPair(2)) * cos(lonLatPair(1));
yPt = cos(lonLatPair(2)) * sin(lonLatPair(1));
zPt = sin(lonLatPair(2));
xyz = [xPt;yPt;zPt];
end