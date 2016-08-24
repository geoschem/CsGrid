function [ Uc, Vc ] = restaggerA2C( Ua, Va )
%RESTAGGERA2C Restagger grids on cubed sphere from A to C (Arakawa grids)

nPerSide = size(Ua,1);
Uc = zeros(nPerSide+1,nPerSide+1,6);
Vc = zeros(nPerSide+1,nPerSide+1,6);


end

