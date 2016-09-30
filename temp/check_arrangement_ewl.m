clear all
close all % close previous figures

% ---load other modules. need to modify the path
CSGridDir = '/n/home08/elundgren/GCHP/tools/CSGrid/';
addpath(genpath(CSGridDir));

% ---grid parameters
NX=48;
Nlon=72;
Nlat=46;

% ---------------------------------
%  generate CS data for testing
%----------------------------------

CSdata = zeros(NX,NX*6);

% a "F"-like graph to check the cube-face orientation. Seems silly.. 
F = zeros(NX,NX);
Fwidth = 2; 
F(fix(NX/3)-Fwidth:fix(NX/3)+Fwidth,fix(NX*1/4):fix(NX*3/4))=1;
F(fix(NX/3):fix(NX*2/3),fix(NX/2)-Fwidth:fix(NX/2)+Fwidth)=1;
F(fix(NX/3):fix(NX*2/3),fix(NX*3/4)-Fwidth:fix(NX*3/4)+Fwidth)=1;

% the back ground value equals to the panel number.
% to check the cube-face numbering
for i=1:6
    CSdata(:,i*NX-NX+1:i*NX)=i*(ones(NX,NX)-F); 
end

figure;
for i=1:6
    subplot(2,3,i)
    surf(CSdata(:,i*NX-NX+1:i*NX)','EdgeColor','None');
    colorbar
    caxis([0,6])
    xlim([1,NX])
    ylim([1,NX])
    view(2); 
    title(sprintf('panel %d',i));
end   
suptitle('original cs')

% ----------------------------------------------------
%  regrid cubed sphere to lat/lon using GMAO tile file
% ----------------------------------------------------

xData_gmao = readTileFile('DE0072xPE0046_CF0048x6C.bin');

% first regrid to LL. Use it as the source data.
[ LLdata, xData_gmao ] = regridConservative( CSdata , xData_gmao );

figure;
surf(LLdata','EdgeColor','None');
colorbar
caxis([0,6])
xlim([1,Nlon])
ylim([1,Nlat])
view(2); 
caxis([0,6])
suptitle('original cs --gmao--> lat/lon')

% ----------------------------------------------------
%  regrid lat/lon back to cubed sphere
% ----------------------------------------------------

% Then regrid back to CS. The result matches the original CS data
[ CSdata_GMAO, xData_gmao ] = regridConservative( LLdata , xData_gmao );

figure;
for i=1:6
    subplot(2,3,i)
    surf(CSdata_GMAO(:,i*NX-NX+1:i*NX)','EdgeColor','None');
    colorbar
    caxis([0,6])
    xlim([1,NX])
    ylim([1,NX])
    view(2); 
    title(sprintf('panel %d',i));
end  
suptitle('lat/lon --gmao--> cs')

% ---------------------------------------------------------
%  read in tempest netcdf mapping file (lat/lon -> cs)
% ---------------------------------------------------------

L2C_file='4x5-to-c48_MAP.nc';
L2C_L_1D = ncread(L2C_file,'col');
L2C_C_1D = ncread(L2C_file,'row');
L2C_S = ncread(L2C_file,'S');

% -----------------------------------------------------------------
%  regrid lat/lon to cubed sphere with Tempest MAP weights (Jaiwei)
% -----------------------------------------------------------------
L2C_L = convert_1Dto2D(L2C_L_1D,Nlon,Nlat);
L2C_C = convert_1Dto2D(L2C_C_1D,NX,NX*6);

yData_j = xData_gmao; %get the same format
yData_j(1).W=L2C_S; % 1=LL and C2L; 2=CS and L2C
yData_j(1).II=L2C_L(:,1);
yData_j(1).JJ=L2C_L(:,2);
yData_j(2).II=L2C_C(:,1);
yData_j(2).JJ=L2C_C(:,2);
yData_j(2).W=L2C_S; % for C2L we need to rearrange this!

[ CSdata_tempest_j, yData_j ] = regridConservative( LLdata , yData_j );

figure;
for i=1:6
    subplot(2,3,i)
    surf(CSdata_tempest_j(:,i*NX-NX+1:i*NX)','EdgeColor','None');
    colorbar
    caxis([0,6])
    xlim([1,NX])
    ylim([1,NX])
    view(2); 
    title(sprintf('panel %d',i));
end   
suptitle('lat/lon --tempest--> cs (Jaiwei)')

% ---------------------------------------------------------
%  Lizzie's Fix #1 - Fix the 1D to 2D mapping
% ---------------------------------------------------------

% Problem: some L2C_L(:,2) and L2C_C(:,2) values are off by 1

% Replace convert_1Dto2D with built-in matlab function 'ind2sub'
numPoints = length(L2C_L_1D);
llarrSize = [Nlon,Nlat];
csarrSize = [NX,NX*6];
L2C_L_l = zeros(numPoints,2);
L2C_C_l = zeros(numPoints,2);
for i=1:numPoints
    [L2C_L_l(i,1), L2C_L_l(i,2)] = ind2sub(llarrSize,L2C_L_1D(i));
    [L2C_C_l(i,1), L2C_C_l(i,2)] = ind2sub(csarrSize,L2C_C_1D(i));
end

yData_l = xData_gmao; %get the same format
yData_l(1).II = L2C_L_l(:,1);
yData_l(1).JJ = L2C_L_l(:,2);
yData_l(2).II = L2C_C_l(:,1);
yData_l(2).JJ = L2C_C_l(:,2);
yData_l(1).W = L2C_S;
yData_l(2).W = L2C_S;

[ CSdata_tempest, yData_l ] = regridConservative( LLdata , yData_l );
figure;
for i=1:6
    subplot(2,3,i)
    surf(CSdata_tempest(:,i*NX-NX+1:i*NX)','EdgeColor','None');
    colorbar
    caxis([0,6])
    xlim([1,NX])
    ylim([1,NX])
    view(2); 
    title(sprintf('panel %d',i));
end   
suptitle('lat/lon --tempest--> cs (Lizzie Fix 1 - 1D to 2D fix)')

% --------------------------------------------------------------
%  Lizzie's Fix #2: correct longitude shift
% --------------------------------------------------------------

llII = yData_l(1).II;
llII_new = zeros(length(llII),1);
ii_offset = -2;
for i = 1:length(llII);
    if ii_offset > 0 &  llII(i) > (Nlon-ii_offset);
        llII_new(i) = llII(i) + ii_offset - Nlon;
    elseif ii_offset < 0 & llII(i) < abs(ii_offset)+1; 
        llII_new(i) = llII(i) + ii_offset + Nlon;
    else
        llII_new(i) = llII(i) + ii_offset;
    end
end

% Make new yData from last generated yData, and update with this fix
yData_l2 = yData_l;
yData_l2(1).II=llII_new;

[ CSdata_tempest, yData_l2 ] = regridConservative( LLdata , yData_l2 );
figure;
for i=1:6
    subplot(2,3,i)
    surf(CSdata_tempest(:,i*NX-NX+1:i*NX)','EdgeColor','None');
    colorbar
    caxis([0,6])
    xlim([1,NX])
    ylim([1,NX])
    view(2); 
    title(sprintf('panel %d',i));
end   
suptitle('lat/lon --tempest--> cs (Lizzie Fix 2 - longitude offset fix)')

% --------------------------------------------------------------
%  Lizzie's Fix #3: remap CS JJ to map to correct faces
% --------------------------------------------------------------
csJJ = yData_l2(2).JJ;
csJJ_new = zeros(length(csJJ),1);

% loop over cubed sphere faces
for f=1:6;

   % define new blocks for faces (f is orig, a is after)
   if f==1; a=4; end; % change 1:48 to 145:192
   if f==2; a=5; end; % etc
   if f==3; a=1; end;
   if f==4; a=2; end;
   if f==5; a=6; end;
   if f==6; a=3; end;

   % define the face block start and end indexes (e.g. 1 and 48 for f=1)
   face_start_ind_orig = NX*(f-1)+1;
   face_end_ind_orig = NX*f;
   face_start_ind_new = NX*(a-1)+1; % e.g. 145 for first original face (f==1)

   % get JJ indexes corresponding to this face
   this_face_ind = find(csJJ >= face_start_ind_orig & ...
			csJJ <= face_end_ind_orig);

   % Loop over all JJ indexes for this face		       
   for j = 1:length(this_face_ind);

       % determine the offset of this JJ value from original face block start
       offset_face_JJ = csJJ(this_face_ind(j)) - face_start_ind_orig;

       % set a new JJ value that is the same offset from new face block start
       csJJ_new(this_face_ind(j)) = face_start_ind_new + offset_face_JJ;
   end 
end

% Make new yData from last generated yData, and update
yData_l3 = yData_l2;
yData_l3(2).JJ=csJJ_new;

[ CSdata_tempest, yData_l3 ] = regridConservative( LLdata , yData_l3 );
figure;
for i=1:6
    subplot(2,3,i)
    surf(CSdata_tempest(:,i*NX-NX+1:i*NX)','EdgeColor','None');
    colorbar
    caxis([0,6])
    xlim([1,NX])
    ylim([1,NX])
    view(2); 
    title(sprintf('panel %d',i));
end   
suptitle('lat/lon --tempest--> cs (Lizzie Fix 3 - map to new faces)')

% --------------------------------------------------------------
%  Lizzie's Fix #4: rotation issues
% --------------------------------------------------------------
csII = yData_l3(2).II;
csJJ = yData_l3(2).JJ;
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
                    (face_IImax_ind - csII_new(this_face_ind(j)) );
		                          
       end	  


   end 
end

% Make new yData from last generated yData, and update
yData_l4 = yData_l3;
yData_l4(2).II=csII_new;
yData_l4(2).JJ=csJJ_new;

[ CSdata_tempest, yData_l4 ] = regridConservative( LLdata , yData_l4 );
figure;
for i=1:6
    subplot(2,3,i)
    surf(CSdata_tempest(:,i*NX-NX+1:i*NX)','EdgeColor','None');
    colorbar
    caxis([0,6])
    xlim([1,NX])
    ylim([1,NX])
    view(2); 
    title(sprintf('panel %d',i));
end   
suptitle('lat/lon --tempest--> cs (Lizzie Fix 4 - rotation issues)')


% ---------------------------------------------------------
%  Write stuff to files for inspection
% ---------------------------------------------------------

gmao_cs_II = xData_gmao(2).II;
gmao_cs_JJ = xData_gmao(2).JJ;
gmao_C_2D = [gmao_cs_II gmao_cs_JJ];
dlmwrite('gmao_C_2D.csv',gmao_C_2D);

gmao_ll_II = xData_gmao(1).II;
gmao_ll_JJ = xData_gmao(1).JJ;
gmao_L_2D = [gmao_ll_II gmao_ll_JJ];
dlmwrite('gmao_L_2D.csv',gmao_L_2D);

% ---------------------------------------------------------
%  Preliminary C2L code - not yet used
% ---------------------------------------------------------

C2L_file='c48-to-4x5_MAP.nc';
C2L_C_1D = ncread(C2L_file,'col');
C2L_L_1D = ncread(C2L_file,'row');
C2L_S = ncread(C2L_file,'S');

C2L_L = convert_1Dto2D(C2L_L_1D,Nlon,Nlat);
C2L_C = convert_1Dto2D(C2L_C_1D,NX,NX*6);


