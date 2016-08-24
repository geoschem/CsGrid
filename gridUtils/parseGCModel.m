function [ vertGridRes ] = parseGCModel( vertGridRes )
%PARSEGCMODEL Convert a GEOS-Chem BPCH model string to a vertical res.
%   GEOS-Chem BPCH files describe their vertical resolution using a single
%   string, which must be parsed. This function will convert that string
%   into vertGridRes, which is a resolution string which can be processed
%   by parseGridVert.

OKModels = {'GEOSFP','MERRA','GEOS5',...
    'GEOSFP_47L','MERRA_47L','GEOS5_47L',...
    'GISSM23','GISSF23','GISSX23',...
    'GISSM40','GISSF40','GISSX40',...
    'GISSM102','GISSF102','GISSX102'};
    
if nargin < 1
    % Output a list of OK resolutions
    vertGridRes = [];
    fprintf('Known GC model strings:\n');
    for iModel = 1:length(OKModels)
        fprintf(' => %s\n',OKModels{iModel});
    end
    return;
elseif strcmpi(vertGridRes,'listmodels')
    % Return a list of OK resolutions
    vertGridRes = OKModels;
    return;
else
    % Allow for special cases with enhanced vertical grids
    if any(strncmpi(vertGridRes,{'GEOSFP','MERRA2','MERRA','GEOS5'},5))
        resSplit = strsplit(vertGridRes,'_');
        if length(resSplit) == 1
            % No custom size
            vertGridRes = 'GEOS5';
        else
            nLev = str2double(resSplit{2}(1:(end-1)));
            if nLev == 47
                % Reduced vertical grid
                vertGridRes = 'GEOS5RV';
            elseif nLev == 72
                % Standard vertical grid
                vertGridRes = 'GEOS5';
            else
                % Enhanced vertical grid
                vertGridRes = sprintf('GEOS5L%i',nLev);
            end
        end
    else
        switch upper(vertGridRes)
            case {'EXTERNAL'}
                % Nothing to be done
                warning('parseGCModel:externalRes','External grids should be specified independently');
                %{
            case {'GEOSFP','MERRA','GEOS5','MERRA2'}
                % All have the same vertical grid
                vertGridRes = 'GEOS5';
            case {'GEOSFP_47L','MERRA_47L','GEOS5_47L','MERRA2_47L'}
                % Same reduced vertical grid
                vertGridRes = 'GEOS5RV';
                %}
            case {'GISSM23','GISSF23','GISSX23'}
                % GISS (ModelE) 23-layer
                vertGridRes = 'GISS23';
            case {'GISSM40','GISSF40','GISSX40'}
                % GISS (ModelE) 40-layer
                vertGridRes = 'GISS40';
            case {'GISSM102','GISSF102','GISSX102'}
                % GISS (ModelE) 102-layer
                vertGridRes = 'GISS102';
            case {'GEOS4','GCAP','GEOS5_30L'}
                warning('parseGCModel:obsoleteRes','Model ''%s'' is obsolete and will not be correctly processed',vertGridRes);
            otherwise
                error('parseGCModel:badRes','Model ''%s'' not recognized',vertGridRes);
        end
    end
end

end

