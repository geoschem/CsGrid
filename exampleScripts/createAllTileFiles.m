inputDir = '/acenet/shared/clee/offline_MAP';
outputDir = 'temp';

file_list = dir(inputDir);
for ii = 1:length(file_list)
    tempestPath = file_list(ii).folder;
    tempestFile = file_list(ii).name;
    tokens = regexp(tempestFile, 'lon(\d+)_lat(\d+)-to-c(\d+)_MAP_(D[CE])(P[CE])_([a-zA-Z\d]+)_?([a-zA-Z\d]+)?.nc', 'tokens');
    if ~isempty(tokens)
        tokens = tokens{1};
        Nlon = str2num(tokens{1});
        Nlat = str2num(tokens{2});
        NX = str2num(tokens{3});
        dateline = tokens{4};
        polar = tokens{5};
        % for now I'm not sure what to do with token 6 since it's always
        % GMAO. Just ignoring.
        
        xData_temp = readTempest( [ tempestPath '/' tempestFile ], Nlon, Nlat, NX, dateline, polar);
        % if there's more information after the GMAO then we are looking
        % at a regional file and should give it the UU prefixes
        if (isempty(tokens{7}))
            TileFileName = sprintf('%s%04dx%s%04d_CF%04dx6C.bin', ...
                dateline, Nlon, polar, Nlat, NX);
        else
            xData_temp
            TileFileName = sprintf('UU%04dxUU%04d_CF%04dx6C.bin', ...
                Nlon, Nlat, NX);
            % The grid names in the tilefile must match the filename
            xData_temp(1).Name = sprintf('UU%04dxUU%04d', Nlon, Nlat);
        end
        
        writeTileFile([outputDir '/' TileFileName], xData_temp);
        xData_new = displayTileFile( [outputDir '/' TileFileName] );
    end
end
