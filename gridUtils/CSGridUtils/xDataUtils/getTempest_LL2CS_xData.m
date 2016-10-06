function [ xData_temp ] = getTempest_LL2CS_xData( tempfile, LLname, CSname, ...
					        Nlat, Nlon, NX )
% GETTEMPESTXDATA Read in Tempest lat/lon to cubed sphere netcdf file and
%                 create an xData structure array that can produce a tilefile
%
%    Input:  (1) tempfile: Tempest netcdf file (lat/lon -> CS only!)   
%
%    Output: (1) xDataOut: xData struct
%
%    Use Case:
%
% Lizzie Lundgren, 10/6/16

%  read in tempest netcdf mapping file (lat/lon -> cs)
L2C_L_1D = ncread( tempfile, 'col');
L2C_C_1D = ncread( tempfile, 'row');
L2C_S    = ncread( tempfile, 'S'  );

llSize = [Nlon,Nlat];
csSize = [NX,NX*6];
numPoints = length(L2C_L_1D);

% Initialize tempest xData structure, using GMAO values except for mapping info
xData_temp = struct('Name', {LLname; CSname}, ...
	            'NX',   {Nlon; NX},       ...
		    'NY',   {Nlat; NX*6},     ...
		    'II',   cell(2,1),        ...
		    'JJ',   cell(2,1),        ...
                    'W',    cell(2,1));

% Set mapping variables
for i=1:numPoints
    [xData_temp(1).II(i), xData_temp(1).JJ(i)] = ind2sub(llSize,L2C_L_1D(i));
    [xData_temp(2).II(i), xData_temp(2).JJ(i)] = ind2sub(csSize,L2C_C_1D(i));
end
xData_temp(1).W = L2C_S;
xData_temp(2).W = L2C_S;

% Eventually add the xData corrections to this function!
% More testing is needed first.