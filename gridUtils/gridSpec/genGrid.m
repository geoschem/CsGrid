function [ lonEdge, latEdge ] = genGrid( lonStride,latStride,halfPolar,centre180,offset180 )
%GENGRID Generates lat/lon edges for a grid specification
%   lonStride:      Degrees longitude per cell
%   latStride:      Degrees latitude per cell
%   halfPolar:      Boolean. Are the N/S polar cells half-cells?
%                   Can also be an integer (e.g. 2 half-polar cells)
%   centre180:      Boolean. Is -180 an edge or a grid center?
%   offset180:      Longitude offset value (defaults to zero)

if nargin < 5
    offset180 = 0;
end

if islogical(halfPolar)
    if halfPolar
        halfPolar = 1;
    else
        halfPolar = 0;
    end
end

halfPolar = cast(halfPolar,'like',latStride);

% Latitude
if halfPolar
    latMinP = -90.0 + (halfPolar.*latStride/2.0);
    nLatInt = round((180.0/latStride))-halfPolar;
    nLatTotal = nLatInt + (2.*halfPolar);
    latEdge = [linspace(-90,latMinP,halfPolar+1),...
               linspace(latMinP+latStride,-(latMinP+latStride),nLatInt-1),...
               linspace(-latMinP,90,halfPolar+1)];
else
    nLat = (round(180.0/latStride));
    latEdge = linspace(-90,90,nLat+1);
end
%nLat = length(latEdge) - 1;

% Longitude
nLon = (round(360.0/lonStride));
if centre180
    lonStart = -180.0 - (lonStride/2.0);
else
    lonStart = -180;
end
lonEdge = linspace(lonStart,lonStart+(nLon*lonStride),nLon+1);
lonEdge = lonEdge + offset180;

end

