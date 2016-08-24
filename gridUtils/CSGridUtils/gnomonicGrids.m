function [ lambdaRad, thetaRad ] = gnomonicGrids( gridType, nPerSide )
%GNOMONICGRIDS Adapted from gnomonic_grids in fvcore_grid_utils.F90
% Calculates the locations of the cubed-sphere grid edges in spherical
% co-ordinates (lambda -> longitude, theta -> latitude). This routine only 
% calculates one face of the grid - use mirrorGrids to get the others.

switch gridType
    case 0
        [lambdaRad,thetaRad] = gnomonicED(nPerSide);
    case 1
        %[lonCtr,latCtr] = gnomonicDist(nPerSide);
        error('gnomonicGrids:distNotImplemented','Grid type ''dist'' not implemented');
    case 2
        %[lonCtr,latCtr] = gnomonicAngl(nPerSide);
        error('gnomonicGrids:anglNotImplemented','Grid type ''angl'' not implemented');
    otherwise
        error('gnomonicGrids:unknownType','Grid type ''%i'' not recognized',gridType);
end

if gridType < 3
    [lambdaRad,thetaRad] = symmED(nPerSide,lambdaRad,thetaRad);
    lambdaRad(:,:) = lambdaRad(:,:) - pi;
end

end

function [lambdaRad,thetaRad] = gnomonicED(nPerSide)

persistent invRt3 aSInvRt3;
if isempty(invRt3)
    invRt3 = 1.0/sqrt(3.0);
    aSInvRt3 = asin(invRt3);
end

% Ranges:
% lambda = [0.75*pi, 1.25*pi]
% theta = [-alpha, alpha]
yDelta = 2.*aSInvRt3 / nPerSide;
nP1 = nPerSide + 1;

lambdaRad = zeros(nP1,nP1);
thetaRad = zeros(nP1,nP1);

% Define East-West edges:
for j=1:nP1
    lambdaRad(1,  j)   	= 0.75*pi;                      % West edge
    lambdaRad(nP1,j)    = 1.25*pi;                      % East edge
    thetaRad(1,   j)    = -aSInvRt3 + yDelta*(j-1);     % West edge
    thetaRad(nP1, j)    = thetaRad(1,j);                % East edge
end

% Get North-South edges by symmetry
for i=2:nPerSide
    [lambdaRad(i,1),thetaRad(i,1)] = LLMirror(lambdaRad(1,1),thetaRad(1,1),...
        lambdaRad(nP1,nP1),thetaRad(nP1,nP1),lambdaRad(1,i),thetaRad(1,i));
    lambdaRad(i,nP1) =  lambdaRad(i,1);
    thetaRad(i,nP1) = -thetaRad(i,1);
end

% Set up X-Y-Z co-ordinates
pp = zeros(3,nPerSide+1,nPerSide+1);
% Set 4 corners
pp(1:3,1,1)     = LL2XYZ2(lambdaRad(1,1)        ,thetaRad(1,1));
pp(1:3,nP1,1)	= LL2XYZ2(lambdaRad(nP1,1)      ,thetaRad(nP1,1));
pp(1:3,1,nP1)	= LL2XYZ2(lambdaRad(1,nP1)      ,thetaRad(1,nP1));
pp(1:3,nP1,nP1) = LL2XYZ2(lambdaRad(nP1,nP1)    ,thetaRad(nP1,nP1));

% Map edges on the sphere back to cube:
% Intersections at x=-rsq3
for j=2:nPerSide
    pp(1:3,1,j) = LL2XYZ2(lambdaRad(1,j),thetaRad(1,j));
    pp(2,1,j)   = -pp(2,1,j)*invRt3/pp(1,1,j);
    pp(3,1,j)   = -pp(3,1,j)*invRt3/pp(1,1,j);
end

 for i=2:nPerSide
    pp(1:3,i,1) = LL2XYZ2(lambdaRad(i,1),thetaRad(i,1));
    pp(2,i,1)   = -pp(2,i,1)*invRt3/pp(1,i,1);
    pp(3,i,1)   = -pp(3,i,1)*invRt3/pp(1,i,1);
 end
 
pp(1,:,:) = -invRt3;
for j=2:nP1
    for i=2:nP1
        % Copy y-z face of the cube along j=1
        pp(2,i,j) = pp(2,i,1);
        % Copy along i=1
        pp(3,i,j) = pp(3,1,j);
    end
end

[lambdaRad,thetaRad] = XYZ2LL3(pp);

end

function [lambdaRad,thetaRad] = symmED(nPerSide,lambdaRad,thetaRad)
% Make grid symmetrical to i=im/2+1
nP1 = nPerSide + 1;
for j=2:nP1
    for i=2:nPerSide
       lambdaRad(i,j) = lambdaRad(i,1);
    end
end

for j=1:nP1
    for i=1:(nPerSide/2)
        iSymm = nPerSide + 2 - i;
       avgPt = 0.5*(lambdaRad(i,j) - lambdaRad(iSymm,j));
       lambdaRad(i,j) = avgPt + pi;
       lambdaRad(iSymm,j) = pi - avgPt;
       avgPt = 0.5*(thetaRad(i,j) + thetaRad(iSymm,j));
       thetaRad(i,j) = avgPt;
       thetaRad(iSymm,j) = avgPt;
    end
end

% Make grid symmetrical to j=im/2+1
for j=1:(nPerSide/2)
    jSymm = nPerSide + 2 - j;
    for i=2:nPerSide
       avgPt = 0.5*(lambdaRad(i,j)+lambdaRad(i,jSymm));
       lambdaRad(i, j) = avgPt;
       lambdaRad(i,jSymm) = avgPt;
       avgPt = 0.5*(thetaRad(i,j)-thetaRad(i,jSymm));
       thetaRad(i, j) =  avgPt;
       thetaRad(i,jSymm) = -avgPt;
    end
end

end

function [lonImg,latImg] = LLMirror(lonMir1,latMir1,lonMir2,latMir2,lonRef,latRef)
% Given the "mirror" as defined by (lonRef1, latRef1), (lonRef2, latRef2),
% and center of the sphere, compute the mirror image of (lonRef, latRef)

xyzRef = LL2XYZ2(lonRef,latRef);
xyzMir1 = LL2XYZ2(lonMir1,latMir1);
xyzMir2 = LL2XYZ2(lonMir2,latMir2);
xyzCross = cross(xyzMir1,xyzMir2);
xyzCross = xyzCross ./ sqrt(sum(xyzCross.*xyzCross));
xyzDot = sum(xyzCross.*xyzRef);
xyzImg = xyzRef - (2.*xyzDot.*xyzCross);

% Convert to lat-lon co-ordinates
[lonImg,latImg,~] = XYZ2LL2(xyzImg);

end