function [ C2LWeight, C2LIdx ] = calcCS2LL( nPerSide, gSpec, showDebug )
%CALCCS2LL Calculates weights for CS -> LL bilinear interpolation

error('calcCS2LL:notReady','Routine does not yet produce correct results!');

% Generate the CS grid description
[~,~,lonCtrCSBase,latCtrCSBase,~,xyzCtrCSBase] = calcCSGrid(nPerSide,'units','radians');

% Expand the grid to include edge values
xyzCtrCS = zeros(3,nPerSide+2,nPerSide+2,6);
lonCtrCS = zeros(nPerSide+2,nPerSide+2,6);
latCtrCS = zeros(nPerSide+2,nPerSide+2,6);
for iTile = 1:6
    xyzTile = zeros(3,nPerSide+2,nPerSide+2);
    for iDim = 1:3
        xyzTile(iDim,:,:) = getHalo(squeeze(xyzCtrCSBase(iDim,:,:,:)),iTile,1);
    end
    xyzCtrCS(:,:,:,iTile) = xyzTile;
    lonCtrCS(:,:,iTile) = getHalo(lonCtrCSBase,iTile,1);
    latCtrCS(:,:,iTile) = getHalo(latCtrCSBase,iTile,1);
end

clear xyzCtrCSBase latCtrCSBase lonCtrCSBase;

% Already have the lat-lon description
nLon = gSpec.nLon;
nLat = gSpec.nLat;
%lonEdgeLL = gSpec.lonEdge.*pi./180;
%latEdgeLL = gSpec.latEdge.*pi./180;
lonCtrLL = gSpec.lonMid.*pi./180;
latCtrLL = gSpec.latMid.*pi./180;
xyzCtrLL = zeros(3,nLon,nLat);
for iLon = 1:nLon
    currLon = lonCtrLL(iLon);
    for jLat = 1:nLat
        currLat = latCtrLL(jLat);
        xyzCtrLL(:,iLon,jLat) = LL2XYZ2(currLon,currLat);
    end
end
dLat = gSpec.latStride * pi / 180;
dLon = gSpec.lonStride * pi / 180;

cellFound = false(nLon,nLat);
shortestDist = zeros(6,1);
C2LIdx = zeros(3,nLon,nLat);
iIter = 0;
allDone = false;
tTimer = tic;
tElapsed = 0;
tStep = 0.2;
storeDist = true;
if storeDist
    GCDistGrid = zeros(nLon,nLat,nPerSide+2,nPerSide+2,6);
    GCDistFound = false(size(GCDistGrid));
end
lonMin = lonCtrLL(1);
hPlot = [];
%pData = zeros(nLon,nLat);
%maxIter = max([10,floor(nPerSide/2)]);
maxIter = 10;
while iIter < maxIter && ~allDone
    iIter = iIter + 1;
    dCubeGrid = zeros(nPerSide,nPerSide,6);
    for iFace = 1:6
        xyzCtrCSFace = xyzCtrCS(:,:,:,iFace);
        lonCtrCSFace = lonCtrCS(:,:,iFace);
        latCtrCSFace = latCtrCS(:,:,iFace);
        shDistFace = shortestDist(iFace);
        dCubeFace = dCubeGrid(:,:,iFace);
        for iY = 1:nPerSide
            iYBase = iY + 1;
            for iX = 1:nPerSide
                iXBase = iX + 1;
                if iIter == 1
                    lonCS = lonCtrCSFace(iXBase,iYBase);
                    latCS = latCtrCSFace(iXBase,iYBase);
                    LL1 = [lonCS,latCS];
                    
                    lonCSP1 = nan;
                    latCSP1 = nan;
                    % Haven't dealt properly with corners..
                    iGuess = 1;
                    iX2V = [1,-1,0,1,0,-1];
                    iY2V = [1,-1,1,0,-1,0];
                    while any(isnan([lonCSP1,latCSP1]))
                        % First guess
                        iX2 = iXBase + iX2V(iGuess);
                        iY2 = iYBase + iY2V(iGuess);
                        lonCSP1 = lonCtrCSFace(iX2,iY2);
                        latCSP1 = latCtrCSFace(iX2,iY2);
                        iGuess = iGuess + 1;
                    end
                    LL2 = [lonCSP1,latCSP1];
                    
                    %{
                    lonCSP1 = lonCtrCSFace(iXBase+1,iYBase+1);
                    latCSP1 = latCtrCSFace(iXBase+1,iYBase+1);
                    LL2 = [lonCSP1,latCSP1];
                    if any(isnan(LL2))
                        lonCSP1 = lonCtrCSFace(iXBase-1,iYBase-1);
                        latCSP1 = latCtrCSFace(iXBase-1,iYBase-1);
                        LL2 = [lonCSP1,latCSP1];
                    end
                    %}
                    
                    dTemp = greatCircleDistance(LL1,LL2,1);
                    dCubeFace(iX,iY) = dTemp;
                else
                    dTemp = dCubeFace(iX,iY);
                end
                dCube = iIter .* dTemp;
                jMin = floor(((latCS-dCube+0.5*pi)/dLat)-iIter+1);
                jMax = ceil(((latCS+dCube+0.5*pi)/dLat)+iIter-1);
                jMin = max([1,   jMin]);
                jMax = min([nLat,jMax]);
                if (jMin == 1 || jMax == nLat)
                    iMin = 1;
                    iMax = nLon;
                else
                    % Number of longitude strides from the middle of the
                    % lowest longitude to the center of the cubed sphere
                    nStride = (lonCS - lonMin)/dLon;
                    % Allow for cycling
                    while (nStride < 0)
                        nStride = nStride + nLon;
                    end
                    while nStride > nLon
                        nStride = nStride - nLon;
                    end
                    % Number of strides per cube cell side length
                    cellStride = abs(dCube) / dLon;
                    iMin = floor(nStride-cellStride-iIter+1);
                    iMax =  ceil(nStride+cellStride+iIter-1);
                    %iMin = floor(((lonCS-dCube)/dLon)-iIter+1);
                    %iMax = ceil(((lonCS+dCube)/dLon)+iIter-1);
                    iMin = max([1,   iMin]);
                    iMax = min([nLon,iMax]);
                end
                for jLat = jMin:jMax
                    for iLon = iMin:iMax
                        % Find nearest cubed-sphere cell for this lat/lon
                        if ~cellFound(iLon,jLat)
                            currLonLL = lonCtrLL(iLon);
                            currLatLL = latCtrLL(jLat);
                            shDistFace = 2*pi;
                            miniFind = false;
                            findIdx = nan(3,1);
                            % Scan through local CS centers and find the
                            % closest one to this lat/lon pairing
                            for iYSub = iY:min([nPerSide,iY+1])
                                iYSubBase = iYSub + 1;
                                for iXSub = iX:min([nPerSide,iX+1])
                                    iXSubBase = iXSub + 1;
                                    if ~storeDist || ~GCDistFound(iLon,jLat,iXSubBase,iYSubBase,iFace)
                                        newDist = greatCircleDistance(...
                                            [currLonLL,currLatLL],...
                                            [lonCtrCSFace(iXSubBase,iYSubBase),...
                                            latCtrCSFace(iXSubBase,iYSubBase)],1);
                                        if storeDist
                                            GCDistFound(iLon,jLat,iXSubBase,iYSubBase,iFace) = true;
                                            GCDistGrid(iLon,jLat,iXSubBase,iYSubBase,iFace) = newDist;
                                        end
                                    else
                                        newDist = GCDistGrid(iLon,jLat,iXSubBase,iYSubBase,iFace);
                                    end
                                    if newDist < shDistFace
                                        miniFind = true;
                                        shDistFace = newDist;
                                        findIdx(1) = iXSubBase;
                                        findIdx(2) = iYSubBase;
                                        findIdx(3) = iFace;
                                    end
                                end
                            end
                            if miniFind
                                [miniFind,miniIdx] =...
                                    getClosestIdx(iLon,jLat,findIdx,xyzCtrCSFace,xyzCtrLL);
                            end
                            if miniFind
                                cellFound(iLon,jLat) = miniFind;
                                C2LIdx(:,iLon,jLat) = miniIdx;
                            else
                                C2LIdx(:,iLon,jLat) = [nan,nan,iFace+0.2];
                            end
                        end
                    end
                end
                if showDebug
                    tNew = toc(tTimer);
                    firstStep = isempty(hPlot);
                    allDone = (iIter > 1) && all(all(cellFound));
                    if firstStep || (tNew - tElapsed) > tStep || allDone
                        if ~firstStep
                            tElapsed = tElapsed + tStep;
                        end
                        dataGrid = squeeze(C2LIdx(3,:,:))';
                        cLim = [-0.5,6.5];
                        %cMap = parula(7);
                        cMap = parula(1000);

                        %pData(iLon,jLat) = shDistFace;
                        %cLim = [0,pi];
                        %cMap = parula(1000);
                        %dataGrid = pData';

                        if firstStep
                            [~,~,hPlot] = plotGrid(dataGrid,'layer');
                            set(gca,'clim',cLim);
                            colormap(gca,cMap);
                            colorbar;
                        else
                            set(hPlot,'cdata',dataGrid);
                        end
                        drawnow;
                    end
                end
            end
        end
        shortestDist(iFace) = shDistFace;
        if iIter == 1
            dCubeGrid(:,:,iFace) = dCubeFace;
        end
    end
    allDone = (iIter > 1) && all(all(cellFound));
    nDone = sum(cellFound(:));
    fprintf('Found %i of %i in %0.1f seconds (%6.2f%%)\n',nDone,nLon*nLat,toc(tTimer),nDone*100/(nLon*nLat));
end

% Do we now need an expensive global search?
if ~allDone
    warning('calcCS2LL:globalSearch','%i cells not matched. Performing global search',numel(cellFound)-nDone);
    for iLon = 1:nLon
        for jLat = 1:nLat
            if ~cellFound(iLon,jLat)
                shDistFace = 2*pi;
                findIdx = nan(3,6);
                for iFace = 1:6
                    miniFind = false;
                    for iX = 1:nPerSide
                        iXBase = iX + 1;
                        for iY = 1:nPerSide
                            iYBase = iY + 1;
                            if ~storeDist || ~GCDistFound(iLon,jLat,iXBase,iYBase,iFace)
                                newDist = greatCircleDistance(...
                                    [lonCtrLL(iLon),latCtrLL(jLat)],...
                                    [lonCtrCSFace(iXBase,iYBase),...
                                    latCtrCSFace(iXBase,iYBase)],1);
                                if storeDist
                                    GCDistFound(iLon,jLat,iXBase,iYBase,iFace) = true;
                                    GCDistGrid(iLon,jLat,iXBase,iYBase,iFace) = newDist;
                                end
                            else
                                newDist = GCDistGrid(iLon,jLat,iXBase,iYBase,iFace);
                            end
                            if newDist < shDistFace
                                miniFind = true;
                                shDistFace = newDist;
                                findIdx(1,iFace) = iXSubBase;
                                findIdx(2,iFace) = iYSubBase;
                                findIdx(3,iFace) = iFace;
                            end
                        end
                    end
                    if miniFind
                        shortestDist(iFace) = shDistFace;
                    else
                        shortestDist(iFace) = +Inf;
                    end
                end
                % Go from shortest to longest
                [~,sortOrder] = sort(shortestDist);
                miniFind = false;
                for iSort = 1:6
                    iFace = sortOrder(iSort);
                    if ~miniFind && ~isinf(shortestDist(iFace))
                        [miniFind,miniIdx] = getClosestIdx(...
                            iLon,jLat,findIdx(:,iFace),xyzCtrCS(:,:,:,iFace),xyzCtrLL);
                        if miniFind
                            cellFound(iLon,jLat) = true;
                            C2LIdx(:,iLon,jLat) = miniIdx;
                        end
                    end
                end
            end
        end
    end
end

assert(allDone,'regridCStoLatLon:matchFailed','Could not determine closest index for %i grid cell(s)',numel(cellFound)-nDone);
C2LWeight = zeros(4,nLon,nLat);
nP1 = nPerSide + 1;
%nP2 = nPerSide + 2;
% C2LWeight stores the OFFSET indices (1:nP2)

for iLon = 1:nLon
    for jLat = 1:nLat
        xyzLL = xyzCtrLL(:,iLon,jLat);
        iX = C2LIdx(1,iLon,jLat);
        iY = C2LIdx(2,iLon,jLat);
        iZ = C2LIdx(3,iLon,jLat);
        if (iX == nP1 && iY == nP1)
            D1 = distToSide(xyzCtrCS(:,iX+1,iY,iZ),xyzCtrCS(:,iX,iY+1,iZ),xyzLL);
            D2 = distToSide(xyzCtrCS(:,iX+1,iY,iZ),xyzCtrCS(:,iX,iY  ,iZ),xyzLL);
            D3 = distToSide(xyzCtrCS(:,iX  ,iY,iZ),xyzCtrCS(:,iX,iY+1,iZ),xyzLL);
            
            C2LWeight(1,iLon,jLat) = D1;
            C2LWeight(2,iLon,jLat) = D2;
            C2LWeight(3,iLon,jLat) = 0.0;
            C2LWeight(4,iLon,jLat) = D3;
            
            C2LSum = sum(C2LWeight(:,iLon,jLat));
            C2LWeight(:,iLon,jLat) = C2LWeight(:,iLon,jLat) ./ C2LSum;
            % Bilinear interp
        elseif (iX == 1 && iY == nP1)
            D1 = distToSide(xyzCtrCS(:,iX+1,iY+1,iZ),xyzCtrCS(:,iX+1,iY+1,iZ),xyzLL);
            D2 = distToSide(xyzCtrCS(:,iX+1,iY  ,iZ),xyzCtrCS(:,iX  ,iY  ,iZ),xyzLL);
            D3 = distToSide(xyzCtrCS(:,iX+1,iY+1,iZ),xyzCtrCS(:,iX  ,iY  ,iZ),xyzLL);
            
            C2LWeight(1,iLon,jLat) = D1;
            C2LWeight(2,iLon,jLat) = 0.0;
            C2LWeight(3,iLon,jLat) = D2;
            C2LWeight(4,iLon,jLat) = D3;
            
            C2LSum = sum(C2LWeight(:,iLon,jLat));
            C2LWeight(:,iLon,jLat) = C2LWeight(:,iLon,jLat) ./ C2LSum;
            % Bilinear interp
        elseif (iX == nP1 && iY == 1)
            % Bilinear interp
            D1 = distToSide(xyzCtrCS(:,iX,iY+1,iZ),xyzCtrCS(:,iX+1,iY+1,iZ),xyzLL);
            D2 = distToSide(xyzCtrCS(:,iX,iY  ,iZ),xyzCtrCS(:,iX+1,iY+1,iZ),xyzLL);
            D3 = distToSide(xyzCtrCS(:,iX,iY  ,iZ),xyzCtrCS(:,iX  ,iY+1,iZ),xyzLL);
            
            C2LWeight(1,iLon,jLat) = D1;
            C2LWeight(2,iLon,jLat) = D2;
            C2LWeight(3,iLon,jLat) = D3;
            C2LWeight(4,iLon,jLat) = 0.0;
            
            C2LSum = sum(C2LWeight(:,iLon,jLat));
            C2LWeight(:,iLon,jLat) = C2LWeight(:,iLon,jLat) ./ C2LSum;
        else
            D1 = distToSide(xyzCtrCS(:,iX  ,iY  ,iZ),xyzCtrCS(:,iX  ,iY+1,iZ),xyzLL);
            D2 = distToSide(xyzCtrCS(:,iX  ,iY+1,iZ),xyzCtrCS(:,iX+1,iY+1,iZ),xyzLL);
            D3 = distToSide(xyzCtrCS(:,iX+1,iY+1,iZ),xyzCtrCS(:,iX+1,iY  ,iZ),xyzLL);
            D4 = distToSide(xyzCtrCS(:,iX+1,iY  ,iZ),xyzCtrCS(:,iX  ,iY  ,iZ),xyzLL);
            
            C2LWeight(1,iLon,jLat) = D2 * D3;
            C2LWeight(2,iLon,jLat) = D3 * D4;
            C2LWeight(3,iLon,jLat) = D4 * D1;
            C2LWeight(4,iLon,jLat) = D1 * D2;
            
            C2LSum = sum(C2LWeight(:,iLon,jLat));
            C2LWeight(:,iLon,jLat) = C2LWeight(:,iLon,jLat) ./ C2LSum;
        end
    end
end

end

function [cellFound,C2LIdx] = getClosestIdx(iLon,jLat,findIdx,xyzCtrFace,xyzCtrLL)
xyzCS = xyzCtrFace;
xyzLL = xyzCtrLL(:,iLon,jLat);
cellFound = false;
C2LIdx = nan(3,1);
iXBase = findIdx(1);
iYBase = findIdx(2);
iX = iXBase - 1;
iY = iYBase - 1;
iZ = findIdx(3); % Face index - don't actually use this
nPerSide = size(xyzCS,2) - 2;
angle1 =  sphereAngle(xyzCS(:,iXBase,iYBase),xyzCS(:,iXBase+1,iYBase  ),xyzCS(:,iXBase,iYBase+1));
angle1a = sphereAngle(xyzCS(:,iXBase,iYBase),xyzCS(:,iXBase+1,iYBase  ),xyzLL);
angle1b = sphereAngle(xyzCS(:,iXBase,iYBase),xyzCS(:,iXBase  ,iYBase+1),xyzLL);
test1  = false;
test11 = false;
test2  = false;
test22 = false;
test3  = false;
test33 = false;
test4  = false;
test44 = false;
if max([angle1a,angle1b]) <= angle1
    test1 = true;
    if (iX == nPerSide && iY == nPerSide)
        angle11  = sphereAngle(xyzCS(:,iXBase+1,iYBase),xyzCS(:,iXBase  ,iYBase+1),xyzCS(:,iXBase,iYBase));
        angle11a = sphereAngle(xyzCS(:,iXBase+1,iYBase),xyzCS(:,iXBase  ,iYBase  ),xyzLL);
        angle11b = sphereAngle(xyzCS(:,iXBase+1,iYBase),xyzCS(:,iXBase  ,iYBase+1),xyzLL);
    else
        angle11  = sphereAngle(xyzCS(:,iXBase+1,iYBase+1),xyzCS(:,iXBase  ,iYBase+1),xyzCS(:,iXBase+1,iYBase));
        angle11a = sphereAngle(xyzCS(:,iXBase+1,iYBase+1),xyzCS(:,iXBase+1,iYBase  ),xyzLL);
        angle11b = sphereAngle(xyzCS(:,iXBase+1,iYBase+1),xyzCS(:,iXBase  ,iYBase+1),xyzLL);
    end
    if max([angle11a,angle11b] <= angle11)
        test11 = true;
        cellFound = true;
        C2LIdx = [iXBase;iYBase;iZ];
    end
else
    angle2 = sphereAngle(xyzCS(:,iXBase,iYBase),xyzCS(:,iXBase,iYBase+1),xyzCS(:,iXBase-1,iYBase));
    angle2a = angle1b;
    angle2b = sphereAngle(xyzCS(:,iXBase,iYBase),xyzCS(:,iXBase-1,iYBase),xyzLL);
    if max([angle2a,angle2b]) <= angle2
        test2 = true;
        if (iX == 1 && iY == nPerSide)
            angle22  = sphereAngle(xyzCS(:,iXBase,iYBase+1),xyzCS(:,iXBase  ,iYBase),xyzCS(:,iXBase-1,iYBase));
            angle22a = sphereAngle(xyzCS(:,iXBase,iYBase+1),xyzCS(:,iXBase-1,iYBase),xyzLL);
            angle22b = sphereAngle(xyzCS(:,iXBase,iYBase+1),xyzCS(:,iXBase  ,iYBase),xyzLL);
        else
            angle22  =sphereAngle(xyzCS(:,iXBase-1,iYBase+1),xyzCS(:,iXBase  ,iYBase+1),xyzCS(:,iXBase-1,iYBase));
            angle22a =sphereAngle(xyzCS(:,iXBase-1,iYBase+1),xyzCS(:,iXBase-1,iYBase  ),xyzLL);
            angle22b =sphereAngle(xyzCS(:,iXBase-1,iYBase+1),xyzCS(:,iXBase  ,iYBase+1),xyzLL);
        end
        if max([angle22a,angle22b] <= angle22)
            test22 = true;
            cellFound = true;
            C2LIdx = [iXBase-1;iYBase;iZ];
        end
    else
        angle3  = sphereAngle(xyzCS(:,iXBase,iYBase),xyzCS(:,iXBase-1,iYBase),xyzCS(:,iXBase,iYBase-1));
        angle3a = angle2b;
        angle3b = sphereAngle(xyzCS(:,iXBase,iYBase),xyzCS(:,iXBase,iYBase-1),xyzLL);
        if max([angle3a,angle3b]) <= angle3 && iX > 1 && iY > 1
            test3 = true;
            angle33  = sphereAngle(xyzCS(:,iXBase-1,iYBase-1),xyzCS(:,iXBase  ,iYBase-1),xyzCS(:,iXBase-1,iYBase));
            angle33a = sphereAngle(xyzCS(:,iXBase-1,iYBase-1),xyzCS(:,iXBase-1,iYBase  ),xyzLL);
            angle33b = sphereAngle(xyzCS(:,iXBase-1,iYBase-1),xyzCS(:,iXBase  ,iYBase-1),xyzLL);
            if max([angle33a,angle33b]) <= angle33
                test33 = true;
                cellFound = true;
                C2LIdx = [iXBase-1;iYBase-1;iZ];
            end
        else
            angle4  = sphereAngle(xyzCS(:,iXBase,iYBase),xyzCS(:,iXBase,iYBase-1),xyzCS(:,iXBase+1,iYBase));
            angle4a = angle3b;
            angle4b = sphereAngle(xyzCS(:,iXBase,iYBase),xyzCS(:,iXBase+1,iYBase),xyzLL);
            if max([angle4a,angle4b]) <= angle4
                test4=true;
                if iX == nPerSide && iY == 1
                    angle44  = sphereAngle(xyzCS(:,iXBase+1,iYBase),xyzCS(:,iXBase,iYBase  ),xyzCS(:,iXBase,iYBase-1));
                    angle44a = sphereAngle(xyzCS(:,iXBase+1,iYBase),xyzCS(:,iXBase,iYBase-1),xyzLL);
                    angle44b = sphereAngle(xyzCS(:,iXBase+1,iYBase),xyzCS(:,iXBase,iYBase  ),xyzLL);
                else
                    angle44  = sphereAngle(xyzCS(:,iXBase+1,iYBase-1),xyzCS(:,iXBase+1,iYBase  ),xyzCS(:,iXBase,iYBase-1));
                    angle44a = sphereAngle(xyzCS(:,iXBase+1,iYBase-1),xyzCS(:,iXBase  ,iYBase-1),xyzLL);
                    angle44b = sphereAngle(xyzCS(:,iXBase+1,iYBase-1),xyzCS(:,iXBase+1,iYBase  ),xyzLL);
                end
                if max([angle44a,angle44b]) <= angle44
                    test44=true;
                    cellFound = true;
                    C2LIdx = [iXBase;iYBase-1;iZ];
                end
            end
        end
    end
end
%if iLon == 127 && jLat >= 57 && jLat <= 59
%    if jLat == 57
%        fprintf('T1? | T2? | T3? | T4?\n');
%    end
%    fprintf('%i %i | %i %i | %i %i | %i %i\n', test1,test11,test2,test22,test3,test33,test4,test44);
%    if jLat == 59
%        pause;
%    end
%end
if cellFound
    assert(C2LIdx(1)>0,'regridCStoLL:lowX','Low X: %i',C2LIdx(1));
    assert(C2LIdx(1)<=(nPerSide+2),'regridCStoLL:highX','High X: %i',C2LIdx(1));
    assert(C2LIdx(2)>0,'regridCStoLL:lowY','Low Y: %i',C2LIdx(2));
    assert(C2LIdx(2)<=(nPerSide+2),'regridCStoLL:highY','High Y: %i',C2LIdx(2));
end
end

function [D2S] = distToSide(v1,v2,pt)
angle = sphereAngle(v1,v2,pt);
side = greatCircleDistance(v1,pt,1);
D2S = asin(sin(side)*sin(angle));
end