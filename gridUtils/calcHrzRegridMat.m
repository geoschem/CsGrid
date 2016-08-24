function [ xMatLon, xMatLat ] = calcHrzRegridMat( gSpecIn, modelOut, makeSparse )
%CALCHRZREGRIDMAT MATLAB implementation of REGRID_A2A_MOD from GEOS-Chem
%   Original routine by S-J Lin.
%   Assumes a regular grid

gSpecOut = gSpecIn.copy;
if ischar(modelOut)
    [ gSpecOut.lonStride,gSpecOut.latStride,gSpecOut.halfPolar,...
        gSpecOut.center180,~,gSpecOut.lonLim,gSpecOut.latLim] ...
        = parseGridHrz( modelOut );
elseif isa(modelOut,'gridSpec')
    %gSpecOut = modelOut;
    gSpecOut.lonStride = modelOut.lonStride;
    gSpecOut.latStride = modelOut.latStride;
    gSpecOut.halfPolar = modelOut.halfPolar;
    gSpecOut.center180 = modelOut.center180;
    gSpecOut.lonLim = modelOut.lonLim;
    gSpecOut.latLim = modelOut.latLim;
    if modelOut.gridSpecial
        gSpecOut.setSpecial(modelOut.lonEdge,modelOut.latEdge);
    end
else
    error('regridHrz:badOutputGrid','Cannot parse output selection');
end
gSpecOut.regenGrid;

nLonIn = gSpecIn.nLon;
nLatIn = gSpecIn.nLat;

% Check the input size
nLonOut = gSpecOut.nLon;
nLatOut = gSpecOut.nLat;
lonEdgeIn = gSpecIn.lonEdge;
latEdgeIn = gSpecIn.latEdge;
lonEdgeOut = gSpecOut.lonEdge;
latEdgeOut = gSpecOut.latEdge;
sinLatIn = sind(latEdgeIn);
sinLatOut = sind(latEdgeOut);

% E-W regridding first
if nLonIn == nLonOut
    xMatLon = eye(nLonIn);
else
    xMatLon = lonRegrid(nLonIn,lonEdgeIn,nLonOut,lonEdgeOut);
end

if nLatIn == nLatOut
    xMatLat = eye(nLatIn);
else
    xMatLat = latRegrid(nLatIn, sinLatIn, nLatOut, sinLatOut);
end

if makeSparse
    xMatLat = sparse(xMatLat);
    xMatLon = sparse(xMatLon);
end

end

function [xMat] = lonRegrid(nLonIn,lonEdgeIn,nLonOut,lonEdgeOut)
%LONREGRID MATLAB implementation of XMAP
% MATLAB does not support FORTRAN-style negative indexing, hence the rather
% tortured indexing
% This function should take one specification and return a matrix, M. For
% each column vector C_In with [N_In x 1] elements, one can then regrid by
%   [      M     ] * [  C_In  ] = [  C_Out  ]
%   [N_Out x N_In] * [N_In x 1] = [N_Out x 1]

xVec = zeros(3*nLonIn + 2,1);
dxVec = zeros(3*nLonIn + 1,1);
xMat = zeros(nLonOut,nLonIn);

xVec((1:(nLonIn+1))+nLonIn+1) = lonEdgeIn;
dxVec((1:nLonIn)+nLonIn+1) = diff(xVec((1+nLonIn)+(1:(nLonIn+1))));

% Ghosting?

% Western edge
found = false;
iLonLow = 1;
while ~found
    if lonEdgeOut(1) >= xVec(iLonLow+nLonIn+1)
        found = true;
    else
        iLonLow = iLonLow - 1;
        if (iLonLow < -nLonIn)
            error('lonregrid:badGhost','Ghosting check failed');
        else
            xVec(iLonLow+nLonIn+1) = xVec(iLonLow+1+nLonIn+1)-dxVec(nLonIn+iLonLow+nLonIn+1);
            dxVec(iLonLow+nLonIn+1) = dxVec(nLonIn+iLonLow+nLonIn+1);
        end
    end
end

% Eastern edge
found = false;
iLonHigh = nLonIn+1;
while ~found
    if lonEdgeOut(end) <= xVec(iLonHigh+nLonIn+1)
        found = true;
    else
        iLonHigh = iLonHigh + 1;
        if (iLonHigh > 2*nLonIn)
            error('lonregrid:badGhost','Ghosting check failed');
        else
            dxVec(iLonHigh+nLonIn) = dxVec(iLonHigh);
            xVec(iLonHigh+nLonIn+1) = xVec(iLonHigh+nLonIn)+dxVec(iLonHigh+nLonIn);
        end
    end
end

dataVec = zeros(3*nLonIn + 1,1,'int32');
dataVec(nLonIn+1) = nLonIn;
dataVec((1:nLonIn)+nLonIn+1) = 1:nLonIn;
dataVec(2*(nLonIn+1)) = 1;
    
% Ghosting?
% Western edge
if iLonLow <= 0
    for iLonIn = iLonLow:0
        dataVec(iLonIn+nLonIn+1) = nLonIn+iLonIn+nLonIn+1;
    end
end

% Eastern edge
if iLonHigh > (nLonIn + 1)
    for iLonIn = (nLonIn+1):(iLonHigh-1)
        dataVec(iLonIn+nLonIn+1) = dataVec(iLonIn+1);
    end
end

iZero = iLonLow;
for iLonOut = 1:nLonOut % i
    found = false;
    iLonIn = iZero;
    while iLonIn <= (iLonHigh-1) && ~found % m
        % Find eastern edge
        found = lonEdgeOut(iLonOut) >= xVec(iLonIn+nLonIn+1) && lonEdgeOut(iLonOut) <= xVec(nLonIn+1+iLonIn+1);
        if found
            noNorm = (lonEdgeOut(iLonOut+1) <= xVec(iLonIn+1+nLonIn+1));
            if noNorm
                % Entire new grid is within original grid
                %dataOut(iLonOut,iLatIn) = dataVec(iLonIn+nLonIn+1);
                xMat(iLonOut,iLonIn) = 1;
                iZero = iLonIn;
            else
                normVal = 1.0/(lonEdgeOut(iLonOut+1)-lonEdgeOut(iLonOut));
                % Left-most fractional area
                %dataSum = (xVec(iLonIn+1+nLonIn+1)-lonEdgeOut(iLonOut))*dataVec(iLonIn+nLonIn+1);
                subLonInMod = mod(iLonIn-1,nLonIn)+1;
                xMat(iLonOut,subLonInMod) = normVal*(xVec(iLonIn+1+nLonIn+1)-lonEdgeOut(iLonOut));
                
                subLonIn = iLonIn+1;
                atEnd = false;
                while subLonIn <= (iLonHigh-1) && ~atEnd % mm
                    atEnd = ~(lonEdgeOut(iLonOut+1) > xVec(subLonIn+1+nLonIn+1));
                    subLonInMod = mod(subLonIn-1,nLonIn)+1;
                    if ~atEnd
                        % Whole cell
                        %dataSum = dataSum + dxVec(subLonIn+nLonIn+1)*dataVec(nLonIn+1+subLonIn);
                        xMat(iLonOut,subLonInMod) = dxVec(subLonIn+nLonIn+1)*normVal;
                    else
                        % Right-most fractional area
                        %dxTemp = lonEdgeOut(iLonOut+1) - xVec(subLonIn+nLonIn+1);
                        xMat(iLonOut,subLonInMod) = (lonEdgeOut(iLonOut+1) - xVec(subLonIn+nLonIn+1))*normVal;
                        %dataSum = dataSum + dxTemp*dataVec(subLonIn+nLonIn+1);
                        iZero = subLonIn;
                    end
                    subLonIn = subLonIn + 1;
                end
            end
        end
        iLonIn = iLonIn + 1;
    end
end

end

function [xMat] = latRegrid(nLatIn, sinLatIn, nLatOut, sinLatOut)

xMat = zeros(nLatOut,nLatIn);
dyVec = diff(sinLatIn);

jZero = 1;
for jLatOut = 1:nLatOut
    jLatIn = jZero;
    sinSouth = sinLatOut(jLatOut);
    found = false;
    while jLatIn <= nLatIn && ~found
        % Find southern edge
        found = (sinSouth>=sinLatIn(jLatIn)) && sinSouth <= sinLatIn(jLatIn+1);
        if found
            noNorm = (sinLatOut(jLatOut+1) <= sinLatIn(jLatIn+1));
            if noNorm
                % Entirely within old gridcell
                %dataOut(:,jLatOut) = dataIn(:,jLatIn);
                xMat(jLatOut,jLatIn) = 1;
                jZero = jLatIn;
            else
                % South-most fractional area
                %dataSum = (sinLatIn(jLatIn+1)-sinLatOut(jLatOut)).*dataIn(:,jLatIn);
                normVal = 1.0/(sinLatOut(jLatOut+1)-sinLatOut(jLatOut));
                xMat(jLatOut,jLatIn) = (sinLatIn(jLatIn+1)-sinLatOut(jLatOut))*normVal;
                subLatIn = jLatIn+1;
                atEnd = false;
                while subLatIn <= nLatIn && ~atEnd
                    atEnd = ~(sinLatOut(jLatOut+1) > sinLatIn(subLatIn+1));
                    if ~atEnd
                        % Intermediate grid cell
                        %dataSum = dataSum + dyVec(subLatIn).*dataIn(:,subLatIn);
                        xMat(jLatOut,subLatIn) = dyVec(subLatIn).*normVal;
                    else
                        % North-most cell
                        dyTemp = sinLatOut(jLatOut+1) - sinLatIn(subLatIn);
                        %dataSum = dataSum + dyTemp.*dataIn(:,subLatIn);
                        xMat(jLatOut,subLatIn) = dyTemp*normVal;
                        jZero = subLatIn;
                    end
                    subLatIn = subLatIn + 1;
                end
            end
        end
        jLatIn = jLatIn + 1;
    end
end

% Average over the polar bands
%{
dataOut(:,1) = sum(dataOut(:,1))/nLonIn;
dataOut(:,end) = sum(dataOut(:,end))/nLonIn;
%}

end