function [lonCorner,latCorner] = validLLPatch(lonCorner,latCorner,mStruct)
%VALIDLLPATCH Re-order lon/lat points to produce a valid quadrilateral
% This code is adapted from the "random_quad" code from the Mathworks blog

if nargin > 2 && ~isempty(mStruct)
    [x,y] = mfwdtran(mStruct,latCorner,lonCorner);
else
    [x,y] = mfwdtran(latCorner,lonCorner);
end
% Triangulate them
t = delaunayTriangulation(x,y);
% If we found 3 triangles, throw one away
if size(t,1) > 2
    t = triangulation(t.ConnectivityList(1:2,:),x,y);
end
% Return the points in the order of the boundary
bound = freeBoundary(t);
bound = bound(:,1);
lonCorner = lonCorner(bound);
latCorner = latCorner(bound);

end