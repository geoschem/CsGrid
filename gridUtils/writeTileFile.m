function [ status ] = writeTileFile( tileFile_new, xData )
% WRITETILEFILE Write a MAPL regrid tile file 
%    Inputs: 
%        (1) Path of the new unformatted FORTRAN binary tilefile
%        (2) xData struct array
%    Output: none
%
% NOTES: The input xData struct array must be the same format as the
%        the xData struct array created by exploreTileFile.m which
%        reads in a binary tilefile. This includes type and x/y
%        locations even though they are not used by MAPL. These variables
%        may be populated with dummy values.
%
% Lizzie Lundgren, 9/26/16

status = 1; % script in progress

fprintf('\n-----------------------------------------------\n');
fprintf('WRITNG NEW BINARY TILEFILE\n');
fprintf('  %s\n', tileFile_new);

%-------------------------------------
% Create and open binary tile file
%-------------------------------------

% Check if file exists. If it does, ask the user whether to overwrite.
% If answer is no, exit.
if exist(tileFile_new) ~= 0;
  fprintf('WARNING: Binary file passed to writeTileFile() already exists.\n');
  answer = input('Proceed and overwrite file? (y/n)\n','s');
  if strcmp(answer,'y');
    delete(tileFile_new);
  else
    assert( strcmp(answer,'n'), 'Not a valid response.');
    fprintf('Exiting writeTileFile().\n');
    status = 2; % file exists. Do not overwrite.
    return
  end
end

% Open file to write
fID = fopen(tileFile_new, 'w');

%------------------------------
% Define data to write
%------------------------------

% Number of data points (cell-to-cell mapping, e.g. length of xData(1).W)
numPoints_memSize = 4;
numPoints = size(xData(1).W,1);

% Number of grids
numGrids_memSize = 4;
numGrids = 2;

% Grid names as 128 byte character arrays (pad with spaces)
gridName_memSize = 128;
grid1_name = [xData(1).Name'; repmat(' ', 128 - length(xData(1).Name), 1)];
grid2_name = [xData(2).Name'; repmat(' ', 128 - length(xData(2).Name), 1)]; 

% NX and NY
numIndexes_memSize = 4;
grid1_NX = xData(1).NX;
grid1_NY = xData(1).NY;
grid2_NX = xData(2).NX;
grid2_NY = xData(2).NY;

% Data arrays
array_memSize = numPoints * numPoints_memSize;
y_dummy = repmat(-1.0, numPoints, 1);
grid1_II = xData(1).II;
grid1_JJ = xData(1).JJ;
grid1_W = xData(1).W;
grid2_II = xData(2).II;
grid2_JJ = xData(2).JJ;
grid2_W = xData(2).W;

% Mystery fields that can be set to dummy values if not present in xData
% since they not used by MAPL. These xData fields are set if reading in
% a GMAO tilefile using displayTileFile.m.

%  (1) Type
if ( isfield(xData,'Type') && isempty(xData(1).Type) ) ...
    | ~isfield(xData(1),'Type');
   type_dummy = repmat(-1.0, numPoints, 1);
else
   type_dummy = xData(1).Type;
end

% (2) X location
if ( isfield(xData,'xLoc') && isempty(xData(1).Xloc) ) ...
   | ~isfield(xData(1),'xLoc');
   x_dummy = repmat(-1.0, numPoints, 1);
else
   x_dummy = xData(1).Xloc;
end

% (3) Y location
if ( isfield(xData,'yLoc') && isempty(xData(1).Yloc) ) ...
   | ~isfield(xData(1),'yLoc');
   y_dummy = repmat(-1.0, numPoints, 1);
else
   y_dummy = xData(1).Yloc;
end

% (4) Mystery data variable (by default, not written to file)
%     NOTE: this is not read by MAPL so can be left out to save space. 
%     However, it is needed if you are replicating a GMAO tilefile and 
%     want identical files (e.g. matching checksums). 
%if isfield(xData,'mysteryVar') && isempty(xData(1).mysteryVar);
%   mystery_var = repmat(-1.0, numPoints, 1);
%else
%   mystery_var= xData(1).mysteryVar;
%end

%------------------------------------------------
% Call local function to write to binary file
%------------------------------------------------

% Number of points
writeVar(fID, numPoints_memSize, numPoints, 'int32');

% Number of grids
writeVar(fID, numGrids_memSize, numGrids, 'int32');

% Grid1 name
writeVar(fID, gridName_memSize, grid1_name, 'char');

% Grid1 # X indexes
writeVar(fID, numIndexes_memSize, grid1_NX, 'int32');

% Grid1 # Y indexes
writeVar(fID, numIndexes_memSize, grid1_NY, 'int32');

% Grid2 name
writeVar(fID, gridName_memSize, grid2_name, 'char');

% Grid2 # X indexes
writeVar(fID, numIndexes_memSize, grid2_NX, 'int32');

% Grid2 # Y indexes
writeVar(fID, numIndexes_memSize, grid2_NY, 'int32');

% Type (not used)
writeVar(fID, array_memSize, type_dummy, 'single');

% X (not used)
writeVar(fID, array_memSize, x_dummy, 'single');

% Y (not used)
writeVar(fID, array_memSize, y_dummy, 'single');

% Grid1 II
writeVar(fID, array_memSize, grid1_II, 'single');

% Grid1 JJ
writeVar(fID, array_memSize, grid1_JJ, 'single');

% Grid1 W
writeVar(fID, array_memSize, grid1_W, 'single');

% Grid2 II
writeVar(fID, array_memSize, grid2_II, 'single');

% Grid2 JJ
writeVar(fID, array_memSize, grid2_JJ, 'single');

% Grid2 W
writeVar(fID, array_memSize, grid2_W, 'single');

%% Mystery variable (turn off by default. see note above)
%writeVar(fID, array_memSize, mystery_var, 'single');

%--------------------------------------------------
% Clean up
%--------------------------------------------------

fclose(fID);
status = 0; % success
end

%--------------------------------------------------
% Function to write variable to binary file
%--------------------------------------------------

function [] = writeVar( fileID, memSize, var, varType )
  fwrite(fileID, memSize, 'int32');  % write mem size [# bytes]
  fwrite(fileID, var, varType);      % write variable
  fwrite(fileID, memSize, 'int32');  % write mem size again for error checking
end 