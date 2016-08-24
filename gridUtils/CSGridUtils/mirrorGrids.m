function [xGrid6,yGrid6] = mirrorGrids(xGrid1,yGrid1)
%MIRRORGRIDS Mirrors one face of a cubed sphere grid onto all 6
rotRadius = 6370.0e3; % in meters, from MAPL_ConstantsMod
nXPt = size(xGrid1,1);
nYPt = size(xGrid1,2);
xGrid6 = zeros(nXPt,nYPt,6);
yGrid6 = zeros(nXPt,nYPt,6);
xGrid6(:,:,1) = xGrid1;
yGrid6(:,:,1) = yGrid1;
iFace = 1;
for j = 1:ceil(nYPt/2.0)
    for i = 1:ceil(nXPt/2.0)
        x1 = 0.25 * (abs(xGrid1(i,j)) + ...
                     abs(xGrid1(nXPt-(i-1),j)) + ...
                     abs(xGrid1(i,nYPt-(j-1))) + ...
                     abs(xGrid1(nXPt-(i-1),nYPt-(j-1))));
                 
        % Be completely certain...
        x1 = abs(x1);
        xGrid6(i,j,iFace)                    = x1.*sign(xGrid1(i,j));
        xGrid6(nXPt-(i-1),j,iFace)           = x1.*sign(xGrid1(nXPt-(i-1),j));
        xGrid6(i,nYPt-(j-1),iFace)           = x1.*sign(xGrid1(i,nYPt-(j-1)));
        xGrid6(nXPt-(i-1),nYPt-(j-1),iFace)  = x1.*sign(xGrid1(nXPt-(i-1),nYPt-(j-1)));
        % Now the latitudes
        y1 = 0.25 * (abs(yGrid1(i,j)) + ...
                     abs(yGrid1(nXPt-(i-1),j)) + ...
                     abs(yGrid1(i,nYPt-(j-1))) + ...
                     abs(yGrid1(nXPt-(i-1),nYPt-(j-1))));
        % Be completely certain...
        y1 = abs(y1);
        yGrid6(i,j,iFace)                    = y1.*sign(yGrid1(i,j));
        yGrid6(nXPt-(i-1),j,iFace)           = y1.*sign(yGrid1(nXPt-(i-1),j));
        yGrid6(i,nYPt-(j-1),iFace)           = y1.*sign(yGrid1(i,nYPt-(j-1)));
        yGrid6(nXPt-(i-1),nYPt-(j-1),iFace)  = y1.*sign(yGrid1(nXPt-(i-1),nYPt-(j-1)));
        %fprintf('i,j,x1,y1: %i %i %f %f\n',i,j,x1,y1);
        % Force consistency on dateline/greenwhich meridian
        if (mod(nXPt,2) ~= 0)
            if ((i == (1+((nXPt-1)/2.0))))
                xGrid6(i,j,iFace) = 0.0;
                xGrid6(i,nYPt-(j-1),iFace) = 0.0;
            end
        end
    end
end

% Now the other faces...
for iFace = 2:6
    for j = 1:nYPt
        for i = 1:nXPt
            x1 = xGrid6(i,j,1);
            y1 = yGrid6(i,j,1);
            z1 = rotRadius;
            xyzIn = [x1;y1;z1];
            if iFace == 2
                % Rotate about z only
                rotAng = -90.0;
                xyzOut = rot3D(3,xyzIn,rotAng);
            elseif iFace == 3
                % Rotate about z, then x
                rotAng = -90.0;
                xyzTemp = rot3D(3,xyzIn,rotAng);
                rotAng = 90.0;
                xyzOut = rot3D(1,xyzTemp,rotAng);
                
                % Force consistency wrt North Pole and dateline/Greenwich
                % meridian
                % These fixes seem to cause distortion of the grid and
                % leave voids?
                %{
                if (mod(nXPt,2)~=0)
                    if ((i==(1+((nXPt-1)/2.0))) && (i==j))
                        xyzOut(1) = 0.0;
                        xyzOut(2) = pi/2.0;
                    end
                    if ((j==(1+((nYPt-1)/2.0))) && (i < (1+((nXPt-1)/2.0))))
                        xyzOut(1) = 0.0;
                    end
                    if ((j==(1+((nYPt-1)/2.0))) && (i > (1+((nXPt-1)/2.0))))
                        xyzOut(1) = pi;
                    end
                end
                %}
            elseif iFace == 4
                rotAng = -180;
                xyzTemp = rot3D(3,xyzIn,rotAng);
                rotAng = 90;
                xyzOut = rot3D(1,xyzTemp,rotAng);
                if ((mod(nXPt,2)~=0) && (j==(1+((nYPt-1)/2.0))))
                    xyzOut(1) = pi;
                end
            elseif iFace == 5
                rotAng = 90;
                xyzTemp = rot3D(3,xyzIn,rotAng);
                rotAng = 90;
                xyzOut = rot3D(2,xyzTemp,rotAng);
            elseif iFace == 6
                rotAng = 90;
                xyzTemp = rot3D(2,xyzIn,rotAng);
                rotAng = 0;
                xyzOut = rot3D(3,xyzTemp,rotAng);
                % Force consistency wrt South Pole and dateline/Greenwich
                % meridian
                % These fixes seem to cause distortion of the grid and
                % leave voids?
                %{
                if (mod(nXPt,2)~=0)
                    if ((i==(1+((nXPt-1)/2.0))) && (i==j))
                        xyzOut(1) = 0.0;
                        xyzOut(2) = -pi/2.0;
                    end
                    if ((i==(1+((nXPt-1)/2.0))) && (j > (1+((nYPt-1)/2.0))))
                        xyzOut(1) = 0.0;
                    end
                    if ((i==(1+((nXPt-1)/2.0))) && (j < (1+((nYPt-1)/2.0))))
                        xyzOut(1) = pi;
                    end
                end
                %}
            end
            xGrid6(i,j,iFace) = xyzOut(1);
            yGrid6(i,j,iFace) = xyzOut(2);
        end
    end
end
end

function [xyzOut] = rot3D(rotAxis,xyzIn,rotAng)
%ROT3D 3-D rotation of a vector about an origin
% rotAxis:  Dimension around which rotation takes place
% xyzIn:    Original vector
% rotAng:   Angle of rotation (degrees)

% Convert spherical to cartesian co-ordinates
xyzCar = zeros(3,1);
[xyzCar(1),xyzCar(2),xyzCar(3)] = sph2cart(xyzIn(1),xyzIn(2),xyzIn(3));

cosAng = cosd(rotAng);
sinAng = sind(rotAng);
switch rotAxis
    case 1
        xyzCar = [  xyzCar(1);...
                    cosAng*xyzCar(2) + sinAng*xyzCar(3);...
                   -sinAng*xyzCar(2) + cosAng*xyzCar(3)];
    case 2
        xyzCar = [  cosAng*xyzCar(1) - sinAng*xyzCar(3);...
                    xyzCar(2);...
                    sinAng*xyzCar(1) + cosAng*xyzCar(3)];
    case 3
        xyzCar = [  cosAng*xyzCar(1) + sinAng*xyzCar(2);...
                   -sinAng*xyzCar(1) + cosAng*xyzCar(2);...
                    xyzCar(3)];
    otherwise
        error('rot3D:badAxis','Axis must be 1, 2 or 3');
end
% Convert back to spherical co-ordinates
xyzOut = zeros(3,1);
[xyzOut(1),xyzOut(2),xyzOut(3)] = cart2sph(xyzCar(1),xyzCar(2),xyzCar(3));
end