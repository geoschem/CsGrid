function [outAngle] = sphereAngle(e1,e2,e3)
% e1: Mid-point
% e2 and e3 to either side
pVec = ones(3,1);
qVec = ones(3,1);
%assert(numel(e1) == 3,'sphereAngle:badE1','First vector does not have 3 elements');
%assert(numel(e2) == 3,'sphereAngle:badE2','Second vector does not have 3 elements');
%assert(numel(e3) == 3,'sphereAngle:badE3','Third vector does not have 3 elements');
pVec(1) = e1(2)*e2(3) - e1(3)*e2(2);
pVec(2) = e1(3)*e2(1) - e1(1)*e2(3);
pVec(3) = e1(1)*e2(2) - e1(2)*e2(1);

qVec(1) = e1(2)*e3(3) - e1(3)*e3(2);
qVec(2) = e1(3)*e3(1) - e1(1)*e3(3);
qVec(3) = e1(1)*e3(2) - e1(2)*e3(1);
ddd = sum(pVec.*pVec) * sum(qVec.*qVec);
if ddd <= 0
    outAngle = 0;
else
    ddd = sum(pVec.*qVec)./sqrt(ddd);
    if (abs(ddd)>1.0)
        % Same thing
        %outAngle = 2*atan(1.0);
        outAngle = pi/2.0;
    else
        outAngle = acos(ddd);
    end
end

end