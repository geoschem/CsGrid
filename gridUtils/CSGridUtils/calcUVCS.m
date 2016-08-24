function [ uLL, vLL ] = calcUVCS( uCS, vCS, noOffset, CStoLL )
%CALCUVCS Calculate lon/lat components of CS winds
% If noOffset, the cube is assumed to have no offset
% CStoLL:   True    Vector is aligned to CS grid, switch to lon/lat
%           False   Vector is lon/lat aligned, switch to CS grid

assert(islogical(CStoLL) && numel(CStoLL) == 1,...
    'calcUVCS:badDirection','Direction must be true (from CS-aligned to LL-aligned) or false (reverse)');

nPerSide = size(uCS,1);

[lonEdge,latEdge,lonMid,latMid] = calcCSGrid(nPerSide,'offsetCube',~noOffset);

if ismatrix(uCS)
    assert(size(uCS,2) == 6*size(uCS,1),'calcUVCS:bad2DSize','2D inputs must be Nx6N');
    % Force 3D
    newU = zeros(nPerSide,nPerSide,6);
    newV = zeros(nPerSide,nPerSide,6);
    for iFace = 1:6
        iLo = (iFace-1)*nPerSide + 1;
        iHi = iLo + nPerSide - 1;
        newU(:,:,iFace) = uCS(:,iLo:iHi);
        newV(:,:,iFace) = vCS(:,iLo:iHi);
    end
    clear uCS vCS;
    uCS = newU;
    vCS = newV;
    clear newU newV;
elseif ndims(uCS) == 3
     assert(size(uCS,3) == 6 && size(uCS,1) == size(uCS,2),'calcUVCS:bad3DSize','3D inputs must be NxNx6');
else
    error('calcUVCS:badDims','Inputs must have 2 or 3 dimensions');
end

uLL = zeros(size(uCS));
vLL = zeros(size(uCS));
d2r = pi/180;

for iFace = 1:6
    u2D = uCS(:,:,iFace);
    v2D = vCS(:,:,iFace);
    lonE2D = d2r.*lonEdge(:,:,iFace);
    latE2D = d2r.*latEdge(:,:,iFace);
    lonM2D = d2r.*lonMid(:,:,iFace);
    latM2D = d2r.*latMid(:,:,iFace);
    for iLon = 1:nPerSide
        for iLat = 1:nPerSide
            LL = [lonE2D(iLon  ,iLat  ),latE2D(iLon  ,iLat  )];
            UL = [lonE2D(iLon  ,iLat+1),latE2D(iLon  ,iLat+1)];
            UR = [lonE2D(iLon+1,iLat+1),latE2D(iLon+1,iLat+1)];
            LR = [lonE2D(iLon+1,iLat  ),latE2D(iLon+1,iLat  )];
            % Need mid-points
            [u2D(iLon,iLat),v2D(iLon,iLat)] = rotateLLCS(u2D(iLon,iLat),...
                v2D(iLon,iLat),midPtSphere(LL,UL),midPtSphere(LL,LR),...
                midPtSphere(LR,UR),midPtSphere(UL,UR),...
                [lonM2D(iLon,iLat),latM2D(iLon,iLat)],~CStoLL);
        end
    end
    uLL(:,:,iFace) = u2D;
    vLL(:,:,iFace) = v2D;
end

end

