function [ pEdge ] = calcPEdge( pOffset, pFactor, pSurf )
%CALCPEDGE Calculate vector of pressure edges
%   Detailed explanation goes here

%if nargin < 3 || all(pSurf == 0) || isempty(pSurf)
if nargin < 3 || isempty(pSurf)
    pSurf = 1013.25;
end

pEdge = pOffset + (pFactor.*pSurf);

end

