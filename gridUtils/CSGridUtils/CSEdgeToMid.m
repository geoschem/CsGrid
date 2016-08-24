function[lonCtr,latCtr,xyzCtr,xyzEdge] = CSEdgeToMid(lonEdge,latEdge)
%CSEDGETOMID Calculate cubed-sphere grid cell centers from edges (rads)
nXPt = size(lonEdge,1) - 1;
nYPt = size(lonEdge,2) - 1;
nFaces = size(lonEdge,3);

xyzEdge = zeros(3,nXPt+1,nYPt+1,nFaces);
xyzCtr = zeros(3,nXPt,nYPt,nFaces);
lonCtr = zeros(nXPt,nYPt,nFaces);
latCtr = zeros(nXPt,nYPt,nFaces);

for iFace = 1:nFaces
    for iXPt = 1:nXPt
        lastX = iXPt == nXPt;
        for iYPt = 1:nYPt
            lastY = iYPt == nYPt;
            % 4 corners
            latCorner = [latEdge(iXPt,iYPt,iFace),latEdge(iXPt+1,iYPt,iFace),...
                        latEdge(iXPt+1,iYPt+1,iFace),latEdge(iXPt,iYPt+1,iFace)];
            lonCorner = [lonEdge(iXPt,iYPt,iFace),lonEdge(iXPt+1,iYPt,iFace),...
                        lonEdge(iXPt+1,iYPt+1,iFace),lonEdge(iXPt,iYPt+1,iFace)];
            % First convert from lat-lon back to cartesian - far too much
            % redundancy here, make more efficient later
            xyzCorner = LL2XYZ2(lonCorner,latCorner);
            
            % Store edge data...
            xyzEdge(:,iXPt,iYPt,iFace) = xyzCorner(:,1);
            if lastX
                xyzEdge(:,iXPt+1,iYPt,iFace) = xyzCorner(:,2);
            end
            if lastX || lastY
                xyzEdge(:,iXPt+1,iYPt+1,iFace) = xyzCorner(:,3);
            end
            if lastY
                xyzEdge(:,iXPt,iYPt+1,iFace) = xyzCorner(:,4);
            end
            
            eMid = sum(xyzCorner,2);
            eAbs = sqrt(sum(eMid.*eMid));
            if eAbs > 0
                eMid = eMid./eAbs;
            end
            xyzCtr(:,iXPt,iYPt,iFace) = eMid;
            [lonCtr(iXPt,iYPt,iFace),latCtr(iXPt,iYPt,iFace)] = ...
                XYZ2LL2(eMid);
        end
    end
end

end