function [ success ] = createCSRestart( input_file, output_file, output_res)
%createCSRestart Regrid a lat-lon restart file into a cubed-sphere one
%   input_file: name of input lat-lon restart file
%   ouptut_file: name of output CS resart file
%   output_res: cs divisions per side (eg. c24 restat would have
%                output_res = 24)

printAll = true;

ncid = netcdf.open(input_file);

try
    if printAll
        fprintf(1, 'Format: %s\n', netcdf.inqFormat(ncid));
    end
    [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
    
    if printAll
        fprintf(1, 'Dimensions\n\tName\t\tLength\n');
    end
    for  dimid = 0:numdims-1
        [dimname, dimlen] = netcdf.inqDim(ncid, dimid);
        if printAll
            fprintf(1, '\t%s\t\t%d\n', dimname, dimlen);
        end
        if (strcmp('lon', dimname))
            lon = dimlen;
            inLonDimId = dimid;
        elseif (strcmp('lat', dimname))
            lat = dimlen;
            inLatDimId = dimid;
        elseif (strcmp('lev', dimname))
            lev = dimlen;
            inLevDimId = dimid;
        elseif (strcmp('time', dimname))
            time = dimlen;
            inTimeDimId = dimid;
        else
            error('Unexpected dimension in input file.');
        end
        
    end
    clear dimid dimname dimlen
    
    
    % let's make sure there's a tilefile for the size we're doing
    %fname = sprintf('/acenet/shared/ctm/GCHP/TileFiles/DC%04dxPC%04d_CF%04dx6C.bin', lon, lat, output_res);
    fname = sprintf('temp/DC%04dxPC%04d_CF%04dx6C.bin', lon, lat, output_res);
    if (~exist(fname, 'file'))
        error('Can''t find appropriate TileFile for %d x %d to c%d', lon, lat, output_res);
    end
    xData = readTileFile(fname);
    
    % set up the new CS restart file
    outid = netcdf.create(output_file, bitor(netcdf.getConstant('NC_NOCLOBBER'), netcdf.getConstant('NC_64BIT_OFFSET')));
    lonDimId = netcdf.defDim(outid, 'lon', output_res);
    latDimId = netcdf.defDim(outid, 'lat', output_res * 6);
    levDimId = netcdf.defDim(outid, 'lev', lev);
    timeDimId = netcdf.defDim(outid, 'time', time);
    spcDimIds = [lonDimId latDimId levDimId timeDimId];
    
    if printAll
        fprintf(1, '\n\nVariables\n\tName\t\tType\tDims\tAttributes\n');
    end
    outVarId = nan(1, numvars);
    for varid = 0:numvars-1
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid, varid);
        if printAll
            [dimname dimlen] = netcdf.inqDim(ncid, dimids(1));
            dimstr = num2str(dimlen);
            for dimid = 2:length(dimids)
                [dimname dimlen] = netcdf.inqDim(ncid, dimids(dimid));
                dimstr = [dimstr 'x' num2str(dimlen)];
            end
            fprintf(1, '%16s\t%d\t%s\t', varname, xtype, dimstr);
        end
        for dimid = 1:length(dimids)
            dimids(dimid) = netcdf.inqDimID(outid, netcdf.inqDim(ncid, dimids(dimid)));
        end
        outVarId(varid+1) = netcdf.defVar(outid, varname, xtype, dimids);
        for attid = 0:natts-1
            attname = netcdf.inqAttName(ncid, varid, attid);
            netcdf.copyAtt(ncid, varid, attname, outid, outVarId(varid+1));
            [xatype attlen] = netcdf.inqAtt(ncid, varid, attname);
            if printAll
                fprintf(1, '%s\t', attname);
            end
        end
        if (printAll)
            fprintf(1, '\n');
        end
        clear attid attname xatype attlen
    end
    if printAll
        fprintf(1, 'Global Attributs\n\tName\t\tLength\n');
    end
    globalVarId = netcdf.getConstant('NC_GLOBAL');
    for attid = 0:numglobalatts-1
        gattname = netcdf.inqAttName(ncid, globalVarId, attid);
        netcdf.copyAtt(ncid, globalVarId, gattname, outid, globalVarId);
        [gxtype gattlen] = netcdf.inqAtt(ncid, globalVarId, gattname);
        if (printAll)
            fprintf(1, '%s\t%d\t%s\n', gattname, gattlen, netcdf.getAtt(ncid, globalVarId, gattname));
        end
    end
    clear attid gattname gxtype gattlen
    
    netcdf.endDef(outid);
    
    % GCHP appears to expect that the restart files hav their levels top
    % down...That is, lev(end) is the lowest and lev(1) is the highest
    % altitude
    levels = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'lev'));
    if ((levels(end) - levels(1)) < 0)
        flipUD = true;
    else
        flipUD = false;
    end
    
    for varid = 0:numvars-1
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid, varid);
        data = netcdf.getVar(ncid, varid);
        if printAll
            fprintf(1, 'Copying variable %s\n', varname);
        end
        % if this is a species variable, need to regrid it
        if (regexp(varname, '^SPC_.*$'))
            fillValue = netcdf.getAtt(ncid, varid, '_FillValue');
            data(data == fillValue) = nan;
            data = regridConservative(data, xData);
            data(isnan(data)) = fillValue;
            if (flipUD)
                data = data(:,:,end:-1:1);
            end
            dimids = spcDimIds;
        elseif strcmp(varname, 'lat')
            data = 1:(output_res*6);
        elseif strcmp(varname, 'lon')
            data = 1:(output_res);
        elseif strcmp(varname, 'lev')
            if (flipUD)
                data = data(end:-1:1);
            end
        end
        % otherwise, just copy the input data
        netcdf.putVar(outid, outVarId(varid+1), data);
    end
    
    clear varid varname xtype dimids natts
    
catch e
    netcdf.close(ncid);
    if (exist('outid', 'var'))
        if (outid > 0)
            netcdf.close(outid);
        end
    end
    rethrow(e);
end
netcdf.close(outid);
netcdf.close(ncid);