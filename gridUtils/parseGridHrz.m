function [ lonStride,latStride,halfPolar,center180,isNested,lonLim,latLim,OKModels,doTranspose] = parseGridHrz( resString )
%PARSEGRIDHRZ Returns longitude and latitude cell edge vectors
%   Supports the following:
%       GMAO4x5         GMAO 4x5
%       GMAO2x25        GMAO 2x2.5
%       GMAO1x125       GMAO 1x1.25
%       GMAO05x0625     GMAO 0.5x0.625
%       GMAO25x03125    GMAO 0.25x0.3125 (GFPNATIVE)
%       GENERIC1x125    Generic 1x1.25
%       GENERIC1x1      Generic 1x1
%       GMAO1x1         GMAO 1x1
%       GFPNATIVE       GEOS-FP native (1/4 x 5/16)
%       [Nested grids]
%       GPW             GPW 1/24x1/24
%       GRUMP           GRUMP 30"x30" (30/3600 deg square)
%       LANDSCAN        LandScan 30"x30"
%       SQR_X           Square grid with cell size X
%       RECT_IxJ        Rectangular grid with lat cell I by lon cell J
%       AUTO_IxJ        Automatic determination (I longitude cells, J lat.)
%       

lonStride = [];
latStride = [];
halfPolar = false;
center180 = false;
% Unless overwritten...
isNested = false;
lonLim = [-Inf,+Inf];
latLim = [-90,90];
doTranspose = false;
OKModels = {'GISSF','GISSM',...
    'GMAO4x5','GMAO2x25','GMAO1x125','GMAO1x1',...
    'GMAO05x0625','GMAO025x03125',...
    'GENERIC1x125','GENERIC1x1','GFPNATIVE',...
    'GRUMP','GPW','LANDSCAN','LS',...
    '05x0666','025x03125',...
    '05x0666na','05x0666us','05x0666eu','05x0666ch',...
    '025x03125na','025x03125us','025x03125se','025x03125ch'};
% Rejected:
%   '025x03125eu'
%   '05x0666se'

% Strip out punctuation
if nargin > 0 && ~isempty(resString)
    miniString = resString;
    miniString(regexp(miniString,'\_')) = [];
    miniString(regexp(miniString,'\.')) = [];
    miniString(regexp(miniString,'-')) = [];
    stdRes = any(strcmpi(OKModels,miniString));
end
if nargin < 1
    % List available models
    fprintf('Pre-set horizonal grids:\n');
    for iModel = 1:length(OKModels)
        fprintf(' => %s\n',OKModels{iModel});
    end
    fprintf('Generic horizontal grids:\n');
    fprintf(' => SQR_I (square IxI)\n');
    fprintf(' => RECT_IxJ (rectangular IxJ)\n');
    fprintf(' => AUTO_IxJ (automatic I lon. cells x J lat. cells)\n');
    return;
elseif isempty(resString)
    % Just wanted list of models
    return;
elseif length(resString) > 4 && ~stdRes
    switch lower(resString(1:4))
        case 'sqr_'
            squareSize = str2double(resString(5:end));
            lonStride = squareSize;
            latStride = squareSize;
            halfPolar = false;
            center180 = false;
        case 'rect'
            resXxY = resString(6:end);
            xLoc = regexpi(resXxY,'x');
            assert(numel(xLoc)==1,'parseGridHrz:badStr','Bad rectangular grid specified');
            % Typically latitude is specified first
            % e.g. 4x5 = 4 deg lat x 5 deg lon
            latStride = str2double(resXxY(1:(xLoc-1)));
            lonStride = str2double(resXxY((xLoc+1):end));
            halfPolar = false;
            center180 = false;
        case 'auto'
            resXxY = resString(6:end);
            xLoc = regexpi(resXxY,'x');
            assert(numel(xLoc)==1,'parseGridHrz:badStr','Bad rectangular grid specified');
            % Typically longitude is specified first
            % e.g. 72x46 = 72 longitude cells by 46 latitude cells
            nLon = str2double(resXxY(1:(xLoc-1)));
            nLat = str2double(resXxY((xLoc+1):end));
            % Build library of possibilities
            nModels = length(OKModels);
            gridFound = false;
            noMatch = false;
            while ~gridFound && ~noMatch
                iModel = 0;
                while ~gridFound && iModel < nModels
                    % Generate horizontal-only grid specification
                    iModel = iModel + 1;
                    gSpec = genGridSpec(OKModels{iModel},[]);
                    % Ignore nested models for now..?
                    %if ~gSpec.isNested
                        resVec = [gSpec.nLon,gSpec.nLat];
                        gridFound = all(resVec == [nLon,nLat]);
                    %end
                end
                if ~gridFound
                    if ~doTranspose
                        nLonOld = nLon;
                        nLon = nLat;
                        nLat = nLonOld;
                        doTranspose = true;
                    else
                        noMatch = true;
                    end
                end
            end
            if gridFound
                if doTranspose
                    warning('parseGridHrz:transposeAutoMatch','Transpose required to match %ix%i grid to library',nLon,nLat);
                end
            else
                error('parseGridHrz:noAutoMatch','Could not match a grid sized %ix%i to any grid in library',nLon,nLat);
            end
            lonStride = gSpec.lonStride;
            latStride = gSpec.latStride;
            halfPolar = gSpec.halfPolar;
            center180 = gSpec.center180;
            % Should be false, but may change this in future
            isNested = gSpec.isNested;
            lonLim = gSpec.lonLim;
            latLim = gSpec.latLim;
        otherwise
            error('parseGridHrz:resUnknown','Resolution ''%s'' not recognised.',resString);
    end
else
    resString = miniString;
    switch (lower(resString))
        case {'gissx'}
            % NASA GISS extra fine (1x1)
            lonStride = 1;
            latStride = 1;
            halfPolar = false;
            center180 = false;
        case {'gissf'}
            % NASA GISS fine (2x2.5)
            lonStride = 2.5;
            latStride = 2;
            halfPolar = false;
            center180 = false;
        case {'gissm'};
            % NASA GISS medium (4x5)
            lonStride = 5;
            latStride = 4;
            halfPolar = false;
            center180 = false;
        case 'gpw'
            lonStride = 1/24;
            latStride = 1/24;
            halfPolar = false;
            center180 = false;
        case {'grump','ls','landscan'}
            lonStride = 1/120;
            latStride = 1/120;
            halfPolar = false;
            center180 = false;
        case {'gmao4x5'}
            lonStride = 5;
            latStride = 4;
            halfPolar = true;
            center180 = true;
        case {'gmao2x25'}
            lonStride = 2.5;
            latStride = 2;
            halfPolar = true;
            center180 = true;
        case {'gmao1x125'}
            lonStride = 1.25;
            latStride = 1;
            halfPolar = true;
            center180 = true;
        case {'gmao05x0625'}
            lonStride = 0.625;
            latStride = 0.5;
            halfPolar = true;
            center180 = true;
        case {'gmao1x1'}
            lonStride = 1;
            latStride = 1;
            halfPolar = true;
            center180 = true;
        case {'generic1x125'}
            lonStride = 1.25;
            latStride = 1;
            halfPolar = false;
            center180 = false;
        case {'generic1x1'}
            lonStride = 1;
            latStride = 1;
            halfPolar = false;
            center180 = false;
        case {'05x0666'}
            lonStride = 2.0/3.0;
            latStride = 0.5;
            halfPolar = true;
            center180 = true;
        case {'gfpnative','025x03125','gmao025x03125'}
            lonStride = 0.3125;
            latStride = 0.25;
            halfPolar = true;
            center180 = true;

            % Nested grids
        case {'05x0666na','05x0666us'}
            lonStride = 2.0/3.0;
            latStride = 0.5;
            halfPolar = true;
            center180 = true;
            isNested = true;
            lonMin = (-140) - (1/3);
            lonMax = ( -40) + (1/3);
            latMin = 9.75;
            latMax = 70.25;
            lonLim = [lonMin,lonMax];
            latLim = [latMin,latMax];
        case '05x0666ch'
            lonStride = 2.0/3.0;
            latStride = 0.5;
            halfPolar = true;
            center180 = true;
            isNested = true;
            lonMin = 70 - (1/3);
            lonMax = 150 + (1/3);
            latMin = -11.25;
            latMax = 55.25;
            lonLim = [lonMin,lonMax];
            latLim = [latMin,latMax];
        case '05x0666eu'
            lonStride = 2.0/3.0;
            latStride = 0.5;
            halfPolar = true;
            center180 = true;
            isNested = true;
            lonMin = -30 - (1/3);
            lonMax = 50 + (1/3);
            latMin = 29.75;
            latMax = 70.25;
            lonLim = [lonMin,lonMax];
            latLim = [latMin,latMax];
        case {'025x03125na','025x03125us'}
            lonStride = 0.3125;
            latStride = 0.25;
            halfPolar = true;
            center180 = true;
            isNested = true;
            lonMin = -130 - (lonStride/2);
            lonMax = -60 + (lonStride/2);
            latMin = 9.75 - (latStride/2);
            latMax = 60 + (latStride/2);
            lonLim = [lonMin,lonMax];
            latLim = [latMin,latMax];
        case '025x03125ch'
            lonStride = 0.3125;
            latStride = 0.25;
            halfPolar = true;
            center180 = true;
            isNested = true;
            lonMin = 70 - (lonStride/2);
            lonMax = 140 + (lonStride/2);
            latMin = 15 - (latStride/2);
            latMax = 55 + (latStride/2);
            lonLim = [lonMin,lonMax];
            latLim = [latMin,latMax];
        case '025x03125se'
            lonStride = 0.3125;
            latStride = 0.25;
            halfPolar = true;
            center180 = true;
            isNested = true;
            lonMin = 75 - (lonStride/2);
            lonMax = 130 + (lonStride/2);
            latMin = 10 - (latStride/2);
            latMax = 30 + (latStride/2);
            lonLim = [lonMin,lonMax];
            latLim = [latMin,latMax];
        case '025x03125eu'
            error('parseGridHrz:EUNestedBad','0.25x0.3125 nested European grid does not make sense. Consult BY');
        otherwise
            error('parseGridHrz:resUnknown','Resolution ''%s'' not recognised.',resString);
    end
end

end
