function [ lonEdge, latEdge, lonCtr, latCtr, xyzEdge, xyzCtr ] = calcCSGrid( nPerSide, varargin )
%CALCCSGRID Calculates the longitudes and latitudes of a cubed-sphere grid
%   Routines are adapted from the "fv_grid_utils.F90" file in GIGC
%   Assuming an equal-distance edge length
% Returns degrees

% Note: only grid type 0 ("true" equal distance) is supported

inParse = inputParser;
addParameter(inParse,'gridType',0,@(x)(validateattributes(x,{'numeric'},...
    {'nonempty','integer','>=',0,'<',3})));
addParameter(inParse,'units','',@(x)(validateattributes(x,{'char'},{})));
addParameter(inParse,'offsetCube',true,@islogical);
parse(inParse,varargin{:});

% Default to "degrees" output
iDeg = 0;
iRad = 1;
if isempty(inParse.Results.units)
    outUnit = iDeg;
else
    switch lower(inParse.Results.units)
        case {'deg','degs','degree','degrees'}
            outUnit = iDeg;
        case {'rad','rads','radian','radians'}
            outUnit = iRad;
        otherwise
            error('calcCSGrid:badUnit','Output unit should be degrees or radians');
    end
end

gridType = inParse.Results.gridType;
offsetCube = inParse.Results.offsetCube;

% Calculate one grid face
[lambdaRad,thetaRad] = gnomonicGrids(gridType,nPerSide);
% Mirror it to the other faces
[lonEdge,latEdge] = mirrorGrids(lambdaRad,thetaRad);

% Clean up the grid
% Shift the corner away from Japan?
if offsetCube
     lonEdge = lonEdge - (pi/18.0);
end
for iFace = 1:6
    for iLon = 1:(nPerSide+1)
        for iLat = 1:(nPerSide+1)
            newVal = lonEdge(iLon,iLat,iFace);
            if newVal < 0
                newVal = newVal + (2*pi);
            end
            if abs(newVal) < 1.0e-10
                newVal = 0.0;
            end
            lonEdge(iLon,iLat,iFace) = newVal;
            if abs(latEdge(iLon,iLat,iFace)) < 1.0e-10
                latEdge(iLon,iLat,iFace) = 0.0;
            end
        end
    end
end

% Clean up the corners?
% This seems to cause more problems than it solves?
%{
lonEdge(1,:,2) = lonEdge(end,:,1);
lonEdge(1,:,3) = lonEdge(end:-1:1,end,1);
lonEdge(:,end,5) = lonEdge(1,end:-1:1,1);
lonEdge(:,end,6) = lonEdge(:,1,1);
lonEdge(:,1,3) = lonEdge(:,end,2);
lonEdge(:,1,4) = lonEdge(end,end:-1:1,2);

lonEdge(end,:,6) = lonEdge(end:-1:1,1,2);
lonEdge(1,:,4) = lonEdge(end,:,3);
lonEdge(1,:,5) = lonEdge(end:-1:1,end,3);
lonEdge(end,:,3) = lonEdge(1,:,4);
lonEdge(:,1,5) = lonEdge(:,end,4);
lonEdge(:,1,6) = lonEdge(end,end:-1:1,4);
lonEdge(1,:,6) = lonEdge(end,:,5);

latEdge(1,:,2) = latEdge(end,:,1);
latEdge(1,:,3) = latEdge(end:-1:1,end,1);
latEdge(:,end,5) = latEdge(1,end:-1:1,1);
latEdge(:,end,6) = latEdge(:,1,1);
latEdge(:,1,3) = latEdge(:,end,2);
latEdge(:,1,4) = latEdge(end,end:-1:1,2);

latEdge(end,:,6) = latEdge(end:-1:1,1,2);
latEdge(1,:,4) = latEdge(end,:,3);
latEdge(1,:,5) = latEdge(end:-1:1,end,3);
latEdge(end,:,3) = latEdge(1,:,4);
latEdge(:,1,5) = latEdge(:,end,4);
latEdge(:,1,6) = latEdge(end,end:-1:1,4);
latEdge(1,:,6) = latEdge(end,:,5);
%}
[lonCtr,latCtr,xyzCtr,xyzEdge] = CSEdgeToMid(lonEdge,latEdge);

% Convert from radians to degrees if necessary
if outUnit == iDeg
    lonCtr = rad2deg(lonCtr);
    latCtr = rad2deg(latCtr);
    lonEdge = rad2deg(lonEdge);
    latEdge = rad2deg(latEdge);
end

end