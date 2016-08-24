classdef gridSpec < handle
    %GRIDSPEC 3D grid specification
    %   Detailed explanation goes here
    
    properties
        % These properties are sufficient to fully define a global or
        % regional rectangular (lat-lon) horizontal grid with a hybrid
        % pressure system
        lonStride=0;
        latStride=0;
        halfPolar=0;
        center180=false;
        offset180=0;
        pOffset=0;
        pFactor=0;
        pSurf=[];
        pMSL = 1013.25; % hPa
        % Nested grids
        lonLim = [];
        latLim = [];
        % Change this at your peril!
        rPlanet = 6.375e6; % in m
    end
    
    properties (SetAccess = private)
        latEdge=[];
        lonEdge=[];
        % Unusual grid?
        gridSpecial=false;
    end
    
    properties (Dependent = true, SetAccess = private)
        % Vectors
        pEdge;
        zEdge;
        pMid;
        zMid;
        lonMid;
        latMid;
        % Calculated using vector methods
        lonMidTrue;
        latMidTrue;
        % Full grid
        pEdge3D;
        zEdge3D;
        pMid3D;
        zMid3D;

        gridArea;
        areaWeights;
        latWeights;
        isNested;
        nLon;
        nLat;
        nLev;
        % Capabilities
        hrzGrid = false;
        vrtGrid = false;
    end
    
    methods
        function gSpec = gridSpec(lonStride_,latStride_,halfPolar_,...
                center180_,pOffset_,pFactor_,lonLim_,latLim_)
            % Do we have horizontal data?
            gSpec.lonStride = lonStride_;
            gSpec.latStride = latStride_;
            gSpec.halfPolar = halfPolar_;
            gSpec.center180 = center180_;
            gSpec.offset180 = 0;
            if nargin < 8
                latLim_ = [-90 90];
                if nargin < 7
                    lonLim_ = [-Inf,+Inf];
                    if nargin < 6
                        pFactor_ = [];
                        if nargin < 5
                            pOffset_ = [];
                        end
                    end
                end
            end
            gSpec.pOffset   = pOffset_;
            gSpec.pFactor   = pFactor_;
            
            gSpec.lonLim    = lonLim_;
            gSpec.latLim    = latLim_;
            gSpec.pSurf     = gSpec.pMSL;
            % Compute grid properties
            gSpec.regenGrid;
        end
        function setSpecial(gSpec,lonEdge_,latEdge_)
            % Special setup
            gSpec.gridSpecial = true;
            gSpec.halfPolar = false;
            gSpec.center180 = false;
            gSpec.lonEdge = lonEdge_;
            gSpec.latEdge = latEdge_;
        end
        function gS2 = copy(gS1)
            gS2 = gridSpec(gS1.lonStride,gS1.latStride,gS1.halfPolar,...
                gS1.center180,gS1.pOffset,gS1.pFactor,gS1.lonLim,gS1.latLim);
            if gS1.gridSpecial
                gS2.setSpecial(gS1.lonEdge,gS1.latEdge);
            end
            gS2.pSurf = gS1.pSurf;
        end
        function disp(gSpec)
            fprintf('Grid specification container\n');
            if gSpec.hrzGrid
                fprintf(' Horizontal specification:\n');
                if gSpec.isNested
                    lonLim_ = gSpec.lonLim;
                    latLim_ = gSpec.latLim;
                    gridDomain = sprintf('%s-%s, %s-%s',...
                        strtrim(getNEWS(lonLim_(1),'lon')),...
                        strtrim(getNEWS(lonLim_(2),'lon')),...
                        strtrim(getNEWS(latLim_(1),'lat')),...
                        strtrim(getNEWS(latLim_(2),'lat')));
                else
                    gridDomain = 'Global';
                end
                fprintf(' => Domain:        %s\n',gridDomain);
                fprintf(' => Resolution:    %-5.2f x %-5.2f\n',gSpec.latStride,gSpec.lonStride);
                fprintf(' => Cells:         %i x %i\n',gSpec.nLon,gSpec.nLat);
            end
            if gSpec.vrtGrid
                fprintf(' Vertical specification:\n');
                fprintf(' => Model top:     %5.4f hPa\n',gSpec.pEdge(end));
                fprintf(' => Levels:        %i\n',gSpec.nLev);
            end
            if ~(gSpec.vrtGrid || gSpec.hrzGrid)
                fprintf(' -- Empty -- ');
            end
        end
            
        function gSpec = regenGrid(gSpec)
            if gSpec.hrzGrid
                if gSpec.gridSpecial
                    lonEdge_ = gSpec.lonEdge;
                    latEdge_ = gSpec.latEdge;
                else
                    [lonEdge_,latEdge_] = genGrid(gSpec.lonStride,...
                        gSpec.latStride,gSpec.halfPolar,...
                        gSpec.center180,gSpec.offset180);
                end
                if gSpec.isNested
                    % Latitude first
                    if diff(gSpec.latLim) < 180
                        latValid = latEdge_>=gSpec.latLim(1) & latEdge_<=gSpec.latLim(2);
                        latEdge_ = latEdge_(latValid);
                    end
                    if diff(gSpec.lonLim) < 360
                        % Longitude is trickier
                        lonLim_ = gSpec.lonLim;
                        if lonLim_(1) > lonLim_(2)
                            lonLim_(1) = lonLim_(1) - 360;
                        end
                        lonEdge_ = [lonEdge_-360,lonEdge_,lonEdge_+360];
                        lonValid = lonEdge_>=lonLim_(1) & lonEdge_<=lonLim_(2);
                        lonEdge_ = lonEdge_(lonValid);
                    end
                end
                gSpec.lonEdge = lonEdge_;
                gSpec.latEdge = latEdge_;
            end
        end
        
        function vrtGrid = get.vrtGrid(gSpec)
            vrtGrid = ~(isempty(gSpec.pOffset) || isempty(gSpec.pFactor));
        end
        
        function hrzGrid = get.hrzGrid(gSpec)
            hrzGrid = ~isempty(gSpec.lonStride+gSpec.latStride);
        end
        
        function lIdx = lonIdx(gSpec,lonVal)
            % Find the index containing the specified longitude
            lIdx = lonVal;
            tooLow = lonVal <= min(gSpec.lonEdge);
            lonVal(tooLow) = lonVal(tooLow) + 360;
            tooHigh = lonVal > max(gSpec.lonEdge);
            lonVal(tooHigh) = lonVal(tooHigh) - 360;
            for iIdx = 1:length(lonVal)
                iVal = lonVal(iIdx);
                lIdx(iIdx) = find(gSpec.lonEdge >= iVal,1,'first') - 1;
            end
        end
        
        function lIdx = latIdx(gSpec,latVal)
            % Find the index containing the specified latitude
            lIdx = latVal;
            for iIdx = 1:length(latVal)
                iVal = latVal(iIdx);
                lIdx(iIdx) = find(gSpec.latEdge >= iVal,1,'first') - 1;
                % Special case
                if lIdx(iIdx) == 0
                    lIdx(iIdx) = 1;
                end
            end
        end
        
        function lIdx = levIdx(gSpec,zVal)
            % Find the index at a specific height (average)
            lIdx = zVal;
            for iIdx = 1:length(zVal)
                iVal = zVal(iIdx);
                lIdx(iIdx) = find(gSpec.zEdge >= iVal,1,'first') - 1;
                % Special case
                if lIdx(iIdx) == 0
                    lIdx(iIdx) = 1;
                end
            end
        end
        
        function [lonIdx,latIdx] = lonLatIdx(gSpec,lonVal,latVal)
            % Find lon/lat indices for a co-ordinate
            lonIdx = gSpec.lonIdx(lonVal);
            latIdx = gSpec.latIdx(latVal);
        end
        
        function nLon = get.nLon(gSpec)
            if gSpec.hrzGrid
                nLon = length(gSpec.lonEdge) - 1;
            else
                nLon = 0;
            end
        end
        
        function nLat = get.nLat(gSpec)
            if gSpec.hrzGrid
                nLat = length(gSpec.latEdge) - 1;
            else
                nLat = 0;
            end
        end
        
        function nLev = get.nLev(gSpec)
            if gSpec.vrtGrid
                nLev = length(gSpec.pFactor) - 1;
            else
                nLev = 0;
            end
        end
        
        function isNested = get.isNested(gSpec)
            isNested = diff(gSpec.latLim) < 180 || diff(gSpec.lonLim) < 360;
        end
        
        function [dx,dy] = dxdy(gSpec)
            % Calculate grid spacing in meters
            rEarth = gSpec.rPlanet;
            latVec = gSpec.latMid;
            % dy (lat spacing) is easy
            dy = rEarth.*diff(latVec).*pi./180.0;
            % dx (lon spacing) is trickier
            % Assume uniform lon spacing (!)
            dx = zeros(size(latVec));
            dLon = diff(gSpec.lonMid(1:2)).*pi./180.0;
            for iLat = 1:gSpec.nLat
                rTemp = rEarth.*sind(abs(90.0-latVec(iLat)));
                dx(iLat) = rTemp.*dLon;
            end
        end
        
        function gArea = get.gridArea(gSpec)
            % Regenerate grid area
            gArea = calcGridArea(gSpec.lonEdge,gSpec.latEdge);
        end
        
        function aWeights = get.areaWeights(gSpec)
            % Get area weights
            gArea = gSpec.gridArea;
            aWeights = gArea./sum(gArea(:));
        end
        
        function lWeights = get.latWeights(gSpec)
            % Get weights by latitude
            gLat = sum(gSpec.gridArea,1);
            lWeights = gLat./sum(gLat);
        end
        
        function pEdge = get.pEdge( gSpec )
            if gSpec.vrtGrid
                pEdge = calcPEdge(gSpec.pOffset, gSpec.pFactor, gSpec.pMSL);
            else
                pEdge = [];
            end
        end

        function zEdge = get.zEdge( gSpec )
            if gSpec.vrtGrid
                zEdge = atmospalt(gSpec.pEdge.*100).*1e-3;
            else
                zEdge = [];
            end
        end
        
        function zMid = get.zMid ( gSpec )
            if gSpec.vrtGrid
                zEdgeV = gSpec.zEdge;
                zMid = (zEdgeV(1:(end-1)) + zEdgeV(2:end))./2;
            else
                zMid = [];
            end
        end
        
        function pMid = get.pMid ( gSpec )
            if gSpec.vrtGrid
                pEdgeV = gSpec.pEdge;
                pMid = (pEdgeV(1:(end-1)) + pEdgeV(2:end))./2;
            else
                pMid = [];
            end
        end
        
        function lonMid = get.lonMid ( gSpec )
            if gSpec.hrzGrid
                lonEdge_ = gSpec.lonEdge;
                lonMid = (lonEdge_(1:(end-1)) + lonEdge_(2:end))./2;
            else
                lonMid = [];
            end
        end
        
        function latMid = get.latMid ( gSpec )
            if gSpec.hrzGrid
                latEdge_ = gSpec.latEdge;
                latMid = (latEdge_(1:(end-1)) + latEdge_(2:end))./2;
            else
                latMid = [];
            end
        end
        
        function lonMidTrue = get.lonMidTrue ( gSpec )
            if gSpec.hrzGrid
                [lonMidTrue,~] = calcMidTrue( gSpec );
            else
                lonMidTrue = [];
            end
        end
        
        function latMidTrue = get.latMidTrue ( gSpec )
            if gSpec.hrzGrid
                [~,latMidTrue] = calcMidTrue( gSpec );
            else
                latMidTrue = [];
            end
        end
        
        function [lonMidTrue,latMidTrue] = calcMidTrue ( gSpec )
            latEdge_ = gSpec.latEdge;
            lonEdge_ = gSpec.lonEdge;
            lonMidTrue = zeros(gSpec.nLon,gSpec.nLat);
            latMidTrue = zeros(gSpec.nLon,gSpec.nLat);
            for iLon_ = 1:gSpec.nLon
                for iLat_ = 1:gSpec.nLat
                    LLCorner = zeros(2,4);
                    LLCorner(:,1) = [lonEdge_(iLon_  );latEdge_(iLat_  )];
                    LLCorner(:,2) = [lonEdge_(iLon_  );latEdge_(iLat_+1)];
                    LLCorner(:,3) = [lonEdge_(iLon_+1);latEdge_(iLat_+1)];
                    LLCorner(:,4) = [lonEdge_(iLon_+1);latEdge_(iLat_  )];
                    LLMid = cellCenter(LLCorner.*pi./180.0).*180.0./pi;
                    lonMidTrue(iLon_,iLat_) = LLMid(1);
                    latMidTrue(iLon_,iLat_) = LLMid(2);
                end
            end
        end
        
        function pEdge = get.pEdge3D( gSpec )
            if gSpec.vrtGrid
                if isempty(gSpec.pSurf)
                    pEdge = repmat(...
                        reshape(gSpec.pEdge,[1,1,gSpec.nLev+1]),...
                        [gSpec.nLon,gSpec.nLat,1]);
                else
                    if gSpec.hrzGrid
                        nLon = gSpec.nLon;
                        nLat = gSpec.nLat;
                    else
                        nLon = size(gSpec.pSurf,1);
                        nLat = size(gSpec.pSurf,2);
                    end
                    pEdge = zeros(nLon,nLat,gSpec.nLev+1);
                    for iLon = 1:nLon
                        for iLat = 1:nLat
                            pEdge(iLon,iLat,:)=calcPEdge(gSpec.pOffset,...
                                gSpec.pFactor, gSpec.pSurf(iLon,iLat));
                        end
                    end
                end
            else
                pEdge = [];
            end
        end

        function zEdge = get.zEdge3D( gSpec )
            if gSpec.vrtGrid
                zEdge = zeros(nLon,nLat,nLev+1);
                pEdge = gSpec.pEdge3D;
                for iLon = 1:nLon 
                    for iLat = 1:nLat
                        pEVec = squeeze(pEdge(iLon,iLat,:));
                        zEdge(iLon,iLat,:) = atmospalt(pEVec.*100).*1e-3;
                    end
                end
            else
                zEdge = [];
            end
        end
        
        function zMid = get.zMid3D ( gSpec )
            if gSpec.vrtGrid
                zEdge = gSpec.zEdge3D;
                zMid = (zEdge(:,:,1:(end-1)) + zEdge(:,:,2:end))./2;
            else
                zMid = [];
            end
        end
        
        function pMid = get.pMid3D ( gSpec )
            if gSpec.vrtGrid
                pEdge = gSpec.pEdge3D;
                pMid = (pEdge(:,:,1:(end-1)) + pEdge(:,:,2:end))./2;
            else
                pMid = [];
            end
        end
        
        function cellDims = cellMetric(gSpec,iLon,iLat,iLev)
            % Calculate cell dimensions in m
            % Vertical first
            cellDims = zeros(1,3);
            zEdgeV = gSpec.zEdge;
            cellDims(3) = 1000.*(zEdgeV(iLev+1) - zEdgeV(iLev));
            latEdgeV = gSpec.latEdge;
            latEdgeV = latEdgeV([iLat,iLat+1]);
            latArc = abs(diff(latEdgeV));
            lonEdgeV = gSpec.lonEdge;
            lonEdgeV = lonEdgeV([iLon,iLon+1]);
            rEarth = gSpec.rPlanet; % in m
            % Latitude edge length is easy
            cellDims(2) = rEarth.*(latArc.*pi./180);
            % Longitude edge we will approximate
            cellArea = GEOSGridArea(latEdgeV,lonEdgeV);
            cellDims(1) = cellArea./cellDims(2);
        end
    end
    
end

function NEWSStr = getNEWS(degVal,dirStr,outFmt)
    switch lower(dirStr)
        case {'lat','latitude'}
            plusStr = 'N';
            negStr = 'S';
        case {'lon','longitude'}
            plusStr = 'E';
            negStr = 'W';
        otherwise
        error('getNEWS:badDir','Bad direction');
    end
    if nargin < 3
        outFmt = '%5.1f%s';
    end
    if degVal < 0
        NEWSStr = sprintf(outFmt,abs(degVal),plusStr);
    else
        NEWSStr = sprintf(outFmt,abs(degVal),negStr);
    end
end
