function [ xData ] = readTempest( tFile, NX, Nlon, Nlat )
% READTEMPEST Read Tempest netcdf file lat/lon to cubed sphere 
%             regridding information to generate an xData struct array 
%             that can generate a GCHP-compatible tile file when passed
%             to writeTileFile().
%
%    Input:  (1) tFile: Tempest output netcdf file (full path)
%            (2) NX:    cubed sphere side length
%            (3) Nlon:  number longitudes
%            (4) Nlat:  number latitudes
%
%    Output: 1) xData:  struct array with same format created when reading
%                       a GMAO tile file using readTileFile.m
%
% NOTES: There is currently a stretching discrepancy between cubed sphere
%        data sets created with GMAO vs Tempest mapping information.
%
% Lizzie Lundgren and Jaiwei Zhang

%------------------------
% Check that file exists
%-----------------------
assert(fastExist(tFile),'readTileFile:fileNotFound', ...
       'Tempest file does not exist');
fID = fopen(tFile,'rb');
cleanupFn = onCleanup(@()(fclose(fID)));

%---------------------------------
% Create xData struct array
%---------------------------------
xCell = cell(2,1);
xData = struct('Name', xCell, ...
               'NX',   xCell, ...
               'NY',   xCell, ...
               'II',   xCell, ...
               'JJ',   xCell, ...
               'W',    xCell);

%-----------
% Read file
%-----------
fprintf('\nREADING TEMPEST FILE:\n');
fprintf('  %s\n', tFile);
L2C_L_1D = ncread(tFile,'col');
L2C_C_1D = ncread(tFile,'row');
L2C_S    = ncread(tFile,'S');

%-------------------------------------
% Convert 1D indexes to 2D subscripts
%-------------------------------------
L2C_L = zeros(length(L2C_L_1D), 2);
L2C_C = zeros(length(L2C_C_1D),2);
for i=1:numPoints
    [L2C_L(i,1), L2C_L(i,2)] = ind2sub([Nlon,Nlat], L2C_L_1D(i));
    [L2C_C(i,1), L2C_C(i,2)] = ind2sub([NX,NX*6], L2C_C_1D(i));
end

%-------------------------------
% Initialize xData
%-------------------------------
% lat/lon
xData(1).Name = 'dummy';
xData(1).NX   = Nlat;
xData(1).NY   = Nlon;
xData(1).II   = L2C_L(:,1);
xData(1).JJ   = L2C_L(:,2);
xData(1).W    = L2C_S;

% cubed sphere
xData(2).Name = 'dummy';
xData(2).NX   = NX;
xData(2).NY   = NX*6;
xData(2).II   = L2C_C(:,1);
xData(2).JJ   = L2C_C(:,2);
xData(2).W    = L2C_S; % is it okay to do this? need to read in C2L data?

%---------------------------
% Correct longitude shift
%---------------------------
indoffset = -2; % need to update this to be generic pi/18 shift (10 deg)
xData(1).II = shiftIndexes(xData(1).II, indOffset, Nlon);

%------------------------------------
% Swap cubed sphere faces
%------------------------------------
newFaceVec = [4 5 1 2 6 3]; % index is current face, value is final face
xData(2).JJ = swapCSFaces( xData(2).JJ, newFaceVec, NX );

%--------------------------
% Transform faces
%--------------------------

% This commented out code will eventually replace the uncommented code below
% Flip both dimensions of face 6 to rotate by 180 degrees
% faces = 6;
%xData(2).II = flipCSFace( xData(2).II, faces, NX );
%xData(2).JJ = flipCSFace( xData(2).JJ, faces, NX);

% Transpose and then flip JJ for faces 3, 4, and 5
% faces = [3 4 5];
%xData(2) = transposeCSFace( xData(2), faces, NX );
%xData(2).JJ = flipCSFace( xData(2).JJ, faces, NX);

csII = xData(2).II;
csJJ = xData(2).JJ;
csII_new = csII;
csJJ_new = csJJ;

% loop over cubed sphere faces
for f=1:6;

   % define new blocks for faces (f is orig, a is after)
   if f==1; action='nothing'; end; % change 1:48 to 145:192
   if f==2; action='nothing'; end; % etc
   if f==3; action='rotate 90deg ccw'; end;
   if f==4; action='rotate 90deg ccw'; end;
   if f==5; action='rotate 90deg ccw'; end;
   if f==6; action='rotate 180deg ccw'; end;

   if strcmp(action,'nothing')
       continue
   end

   % define the face block start and end indexes (e.g. 1 and 48 for f=1)
   face_JJmin_ind = NX*(f-1)+1;
   face_JJmax_ind = NX*f;
   face_IImin_ind = 1;
   face_IImax_ind = NX;

   % get JJ indexes corresponding to this face
   this_face_ind = find(csJJ >= face_JJmin_ind & ...
			  csJJ <= face_JJmax_ind);

   % Loop over all JJ indexes for this face		       
   for j = 1:length(this_face_ind);

       % flip the IIs and JJs for this face
       if strcmp(action,'rotate 180deg ccw')
           csJJ_new(this_face_ind(j)) = face_JJmin_ind + ...
                    (face_JJmax_ind -  csJJ(this_face_ind(j)) );
           csII_new(this_face_ind(j)) = face_IImin_ind + ...
                    ( face_IImax_ind -  csII(this_face_ind(j)) );
		                          
       end	  

       % Transose and then flip the JJs for this face
       if strcmp(action,'rotate 90deg ccw')
           csJJ(this_face_ind(j));
           csII(this_face_ind(j));
           face_JJmin_ind;
           csII_new(this_face_ind(j)) = csJJ(this_face_ind(j)) ...
                                        - face_JJmin_ind + 1;
           csJJ_new(this_face_ind(j)) = face_JJmin_ind ...
                                        + csII(this_face_ind(j)) - 1;
           csII_new(this_face_ind(j)) = face_IImin_ind + ...
                    (face_IImax_ind -  csII_new(this_face_ind(j)) );
		                          
       end	  
   end 
end

% Update xData
xData(2).II=csII_new;
xData(2).JJ=csJJ_new;