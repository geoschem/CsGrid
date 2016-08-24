function [ CSArea ] = calcCSArea( lonEdge, latEdge, rSphere )
%CALCCSAREA Calculate grid area of each CS cell

nPerSide = size(lonEdge,1) - 1;
if nargin < 3 || isempty(rSphere) || rSphere <= 0
    rSphere = 6.375e6; % Radius of earth in m
end
rSq = rSphere*rSphere;

if max(lonEdge(:)) > 6*pi
    %error('calcCSArea:needRads','Edges must be given in radians');
    % Assume they are in degrees
    lonEdge = lonEdge.*pi./180;
    latEdge = latEdge.*pi./180;
end
mStruct = defaultm('globe');
CSArea = zeros(nPerSide,nPerSide,6);
% Need them to be ordered (middle, side, side)
%validCombo = mod([(1:4)',(2:5)',(0:3)'] - 1,4) + 1;
%validCombo = [1,2,4;2,3,1;3,4,2;4,1,3];
validCombo = [1,2,4;2,3,1;3,2,4;4,1,3];
for iFace = 1:6
    faceLon = lonEdge(:,:,iFace);
    faceLat = latEdge(:,:,iFace);
    for iLon = 1:nPerSide
        for iLat = 1:nPerSide
            latCorner = zeros(4,1);
            lonCorner = zeros(4,1);
            xyzVec = zeros(4,3);
            for iVert = 1:4
                xLon = iLon + (iVert > 2);
                xLat = iLat + (iVert == 1 || iVert == 4);
                lonCorner(iVert) = faceLon(xLon,xLat);
                latCorner(iVert) = faceLat(xLon,xLat);
            end
            % Make it into a valid quadrilateral
            %[lonCorner,latCorner] = validLLPatch(lonCorner,latCorner,mStruct);
            % Convert to xyz
            for iVert = 1:4
                xyzVec(iVert,:) = LL2XYZ2(lonCorner(iVert),latCorner(iVert));
            end
            totAng = 0;
            for iCorner = 1:4
                % Get solid angle for the current combination
                currCombo = validCombo(iCorner,:);
                xyzMini = zeros(3,3);
                for iMini = 1:3
                    xyzMini(iMini,:) = xyzVec(currCombo(iMini),:);
                end
                currAng = sphereAngle(xyzMini(1,:),xyzMini(2,:),xyzMini(3,:));
                totAng = totAng + currAng;
            end
            CSArea(iLon,iLat,iFace) = rSq * (totAng - (2.0*pi));
        end
    end
end
end