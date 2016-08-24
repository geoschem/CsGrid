function [ LLData ] = regridCS2LL( CSData, gSpec )
%REGRIDCS2LL Summary of this function goes here
%   Detailed explanation goes here

% Check dimensions
if ndims(CSData) == 3
    % If necessary, reorder the CS data so that it has the original dimensions
    assert(size(CSData,3) == 6 && size(CSData,1) == size(CSData,2),...
        'regridCStoLatLon:badDataSize3D','3D data must be for a single layer only');
elseif ismatrix(CSData)
    assert(size(CSData,2) == 6*size(CSData,1),...
        'regridCStoLatLon:badDataSize2D','2D data must be for a single layer only');
    CSData = reorderCS(CSData,true);
else
    % Invalid
    error('regridCStoLatLon:badDataSizeND','Data must be Nx6N or NxNx6');
end
nPerSide = size(CSData,1);

% Get regridding weights and indices
[C2LWeight, C2LIdx] = calcCS2LL(nPerSide,gSpec,false);

nLon = gSpec.nLon;
nLat = gSpec.nLat;

haloData = zeros(nPerSide+2,nPerSide+2,6);

for iFace = 1:6
    haloData(:,:,iFace) = getHalo(CSData,iFace,1);
end

LLData = zeros(nLon,nLat,'like',CSData);
for iLon = 1:nLon
    for iLat = 1:nLat
        LLMini = 0;
        WVec = C2LWeight(:,iLon,iLat);
        XVec = C2LIdx(1,iLon,iLat) + [0,0,1,1];
        YVec = C2LIdx(2,iLon,iLat) + [0,1,1,0];
        iZ = C2LIdx(3,iLon,iLat);
        for iPt = 1:4
            iX = XVec(iPt);
            iY = YVec(iPt);
            LLMini = LLMini + WVec(iPt) * haloData(iX,iY,iZ);
        end
        LLData(iLon,iLat) = LLMini;
    end
end

end

