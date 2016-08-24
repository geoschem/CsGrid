function [ gSpec, doTranspose ] = genGridSpec( hrzGrid, vrtGrid )
%GENGRIDSPEC Generates a grid specifier object
%   hrzGrid:    Type "parseGridHrz;" to see options
%   vrtGrid:    Type "parseGridVert;" to see options

% parseGridHrz can handle empty input
[lonStride,latStride,halfPolar,centre180,...
    ~,lonLim,latLim,~,doTranspose] = parseGridHrz( hrzGrid );

if nargin < 2
    vrtGrid = [];
end
[pOffset,pFactor] = parseGridVert(vrtGrid);

gSpec = gridSpec(lonStride,latStride,halfPolar,centre180,pOffset,...
    pFactor,lonLim,latLim);

end

