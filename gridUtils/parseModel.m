function [ hrzResName, vertResName ] = parseModel( modelName, resName )
%PARSEMODEL Parses a model name and returns appropriate data
%   Supports the following:
%       modelName       resNames (possible)
%       GEOS-4/5/FP     4x5, 2x2.5, 1x1.25, 1x1 (RV)
%       GISS,MODELE     F102, F53, F40, F23, M40, M23
%       GRUMP           Ignored
%       GPW             Ignored
%       LANDSCAN        Ignored
%       GENERIC         1x1.25, 1x1

RVLoc = regexpi(resName,'rv');
isRV = ~isempty(RVLoc);
% Strip out vertical resolution tags
for iRVLoc = length(RVLoc):-1:1
    resName(RVLoc(iRVLoc) + [0 1]) = [];
end
% Strip out _ and .
resName(regexp(resName,'\_')) = [];
resName(regexp(resName,'\.')) = [];
badRes = false;
switch (lower(modelName))
    case {'giss','modele'}
        % NASA GISS ModelE
        switch upper(resName(1))
            case 'F'
                % Fine (2x2.5)
                hrzResName = 'GISSF';
            case 'M'
                % Medium (4x5)
                hrzResName = 'GISSM';
            case 'C'
                % Coarse (8x10)
                hrzResName = 'GISSC';
            otherwise
                badRes = true;
        end
        vertResName = sprintf('GISS%s',resName(2:end));
    case 'gpw'
        hrzResName = 'GPW';
        vertResName = [];
    case {'landscan','ls','grump'}
        hrzResName = 'GRUMP';
        vertResName = [];
    case {'geos5','geosfp','merra','merra2'}
        vertResName = 'GEOS5';
        if isRV
            vertResName = 'GEOS5RV';
        end
        OKList = {'4x5','2x25','1x125','1x1','05x0666','025x03125',...
            '05x0666na','05x0666us','05x0666ch','05x0666eu',...
            '025x03125na','025x03125us','025x03125ch','025x03125se',...
            '025x03125eu'};
        badRes = ~any(strcmpi(resName,OKList));
        if ~badRes
            hrzResName = sprintf('GMAO%s',lower(resName));
        end
    case {'generic'}
        vertResName = [];
        OKList = {'1x125','1x1'};
        badRes = ~any(strcmpi(resName,OKList));
        if ~badRes
            hrzResName = sprintf('GENERIC%s',lower(resName));
        end
end
if badRes
    error('parseModel:badRes','Resolution ''%s - %s'' not recognized',modelName,resName);
end

end