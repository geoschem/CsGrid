function [ tileData ] = getHalo( CSData, iTile, nHalo )
%GETHALO Returns the 'haloed' data around a given tile
%   Inputs:
%       CSData:     NxNx6 cubed-sphere data, using the MATLAB arrangement
%       iTile:      The tile to be returned
%       nHalo:      Number of halo cells to return (must be <= floor(N/2))

CSDims = size(CSData);
if length(CSDims) == 3
    assert(CSDims(2) == CSDims(1),'getHalo:nonSquareSide','CS data must have square faces (total size NxNx6)');
    nPerSide = CSDims(1);
else
    error('getHalo:badDataDims','CS data must be 3D ([N x N x 6]) in size');
end

inParse = inputParser;
addRequired(inParse,'iTile',@(x)validateattributes(x,{'numeric'},{'integer','scalar','<=',6,'positive'}));
addOptional(inParse,'nHalo',true,@(x)validateattributes(x,{'numeric'},{'integer','<=',floor(nPerSide/2),'nonnegative'}));
parse(inParse,iTile,nHalo);
if nHalo == 0
    warning('getHalo:zeroHalo','No halo requested; returning original tile');
    tileData = CSData(:,:,iTile);
    return;
end
%iStartD = 1;
iStart  = 1 + nHalo;
iEnd    = nPerSide + nHalo;
iEndD   = nPerSide + (nHalo*2);
tileData = nan(iEndD,iEndD);
tileData(iStart:iEnd,iStart:iEnd) = CSData(:,:,iTile);

% Order:
%   1: Tile before start of D1
%   2: Tile before start of D2
%   3: Tile at end of D1
%   4: Tile at end of D2
tileMap = zeros(nPerSide,nHalo,2,4);
switch iTile
    case 1
        borderTiles = [5,3,2,6];
        % Tilemap counts from the lower left, in the output tile's
        % co-ordinate system
        for iHalo = 1:nHalo
            for iCell = 1:nPerSide
                tileMap(iCell,iHalo,:,1) = [iCell,iHalo];
                tileMap(iCell,iHalo,:,2) = [iHalo,iCell];
                tileMap(iCell,iHalo,:,3) = [iHalo,iCell];
                tileMap(iCell,iHalo,:,4) = [iCell,iHalo];
            end
        end
    case 2
        borderTiles = [1,3,4,6];
        for iHalo = 1:nHalo
            for iCell = 1:nPerSide
                tileMap(iCell,iHalo,:,1) = [nPerSide+1-iHalo,iCell];
                tileMap(iCell,iHalo,:,2) = [iCell,nPerSide+1-iHalo];
                tileMap(iCell,iHalo,:,3) = [iCell,nPerSide+1-iHalo];
                tileMap(iCell,iHalo,:,4) = [nPerSide+1-iHalo,iCell];
            end
        end
    case 3
        borderTiles = [1,5,4,2];
        for iHalo = 1:nHalo
            for iCell = 1:nPerSide
                tileMap(iCell,iHalo,:,1) = [iCell,iHalo];
                tileMap(iCell,iHalo,:,2) = [iHalo,iCell];
                tileMap(iCell,iHalo,:,3) = [iHalo,iCell];
                tileMap(iCell,iHalo,:,4) = [iCell,iHalo];
            end
        end
    case 4
        borderTiles = [3,5,6,2];
        for iHalo = 1:nHalo
            for iCell = 1:nPerSide
                tileMap(iCell,iHalo,:,1) = [nPerSide+1-iHalo,iCell];
                tileMap(iCell,iHalo,:,2) = [iCell,nPerSide+1-iHalo];
                tileMap(iCell,iHalo,:,3) = [iCell,nPerSide+1-iHalo];
                tileMap(iCell,iHalo,:,4) = [nPerSide+1-iHalo,iCell];
            end
        end
    case 5
        borderTiles = [3,1,6,4];
        for iHalo = 1:nHalo
            for iCell = 1:nPerSide
                tileMap(iCell,iHalo,:,1) = [iCell,iHalo];
                tileMap(iCell,iHalo,:,2) = [iHalo,iCell];
                tileMap(iCell,iHalo,:,3) = [iHalo,iCell];
                tileMap(iCell,iHalo,:,4) = [iCell,iHalo];
            end
        end
    case 6
        borderTiles = [5,1,2,4];
        for iHalo = 1:nHalo
            for iCell = 1:nPerSide
                tileMap(iCell,iHalo,:,1) = [nPerSide+1-iHalo,iCell];
                tileMap(iCell,iHalo,:,2) = [iCell,nPerSide+1-iHalo];
                tileMap(iCell,iHalo,:,3) = [iCell,nPerSide+1-iHalo];
                tileMap(iCell,iHalo,:,4) = [nPerSide+1-iHalo,iCell];
            end
        end
end

for iCell = 1:nPerSide
    for iHalo = 1:nHalo
        % 'West' face
        iEdge = 1;
        inTile = borderTiles(iEdge);
        d1 = tileMap(iCell,iHalo,1,iEdge);
        d2 = tileMap(iCell,iHalo,2,iEdge);
        tileData(nHalo + 1 - iHalo,nHalo + iCell) = CSData(d1,d2,inTile);
        % 'South' face
        iEdge = 2;
        inTile = borderTiles(iEdge);
        d1 = tileMap(iCell,iHalo,1,iEdge);
        d2 = tileMap(iCell,iHalo,2,iEdge);
        tileData(nHalo + iCell,nHalo + 1 - iHalo) = CSData(d1,d2,inTile);
        % 'East' face
        iEdge = 3;
        inTile = borderTiles(iEdge);
        d1 = tileMap(iCell,iHalo,1,iEdge);
        d2 = tileMap(iCell,iHalo,2,iEdge);
        tileData(nPerSide + nHalo + iHalo,nHalo + iCell) = CSData(d1,d2,inTile);
        % 'North' face
        iEdge = 4;
        inTile = borderTiles(iEdge);
        d1 = tileMap(iCell,iHalo,1,iEdge);
        d2 = tileMap(iCell,iHalo,2,iEdge);
        tileData(nHalo + iCell,nPerSide + nHalo + iHalo) = CSData(d1,d2,inTile);
    end
end

%{
for iCell = 1:nHalo
    for iHalo = 1:nPerSide
        % 'East' face
        iEdge = 2;
        inTile = borderTiles(iEdge);
        tileData(nHalo + 1 - iCell,nHalo + iHalo) = CSData(tileMap(iHalo,iCell,1,iEdge),tileMap(iHalo,iCell,2,iEdge),inTile);
        % 'West' face
        iEdge = 4;
        inTile = borderTiles(iEdge);
        tileData(nPerSide + nHalo + iCell,nHalo + iHalo) = CSData(tileMap(iHalo,iCell,1,iEdge),tileMap(iHalo,iCell,2,iEdge),inTile);
    end
end
%}

end

