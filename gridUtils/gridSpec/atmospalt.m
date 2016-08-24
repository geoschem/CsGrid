function [z_m] = atmospalt(P_Pa)
%ATMOSPALT Lightweight version of atmospalt from the Aerospace Toolbox
% Returns the COESA estimate of altitude (in m) for a pressure (in Pa)

% Generate interpolation vecor
persistent zInt PInt
if isempty(zInt)
    zInt = -1e3:82e3;
    [~,~,PInt] = atmosisa(zInt);
end

%disp(P_Pa)
z_m = interp1(PInt,zInt,P_Pa);

end
