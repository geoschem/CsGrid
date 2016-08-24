function [lonVec,latVec,xyzArr] = XYZ2LL2(xyzArr)
%XYZ2LL2 Convert cartesian (XYZ) to lat-lon co-ordinates
% Also normalizes xyzArr (put on unit sphere)
% Always assume D1 is the target dimension
assert(size(xyzArr,1)==3,'XYZ2LL2:badDims','First dimension must have size 3');
assert(ndims(xyzArr)<=2,'XYZLL2:excessDims','Array must be 2D at most');
nPoints = size(xyzArr,2);
lonVec = zeros(nPoints,1);
latVec = zeros(nPoints,1);
for iPt = 1:nPoints
    xyzVec = xyzArr(:,iPt);
    vecLen = sqrt(sum(xyzVec.*xyzVec));
    xyzVec = xyzVec./vecLen;
    if ((abs(xyzVec(1))+abs(xyzVec(2))) < 1.0e-20)
        lonPt = 0;
    else
        lonPt = atan2(xyzVec(2),xyzVec(1));
    end
    if (lonPt < 0.0)
        lonPt = lonPt + (2.0*pi);
    end
    latPt = asin(xyzVec(3));
    lonVec(iPt) = lonPt;
    latVec(iPt) = latPt;
    % Normalize
    xyzArr(:,iPt) = xyzVec;
end
end