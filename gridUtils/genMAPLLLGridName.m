function gName = genMAPLLLGridName(XYVec,center180,halfPolar)
if nargin == 1
    gName = sprintf('UU%04ixUU%04i',XYVec(1),XYVec(2));
else
    if halfPolar
        poleChar = 'C';
    else
        poleChar = 'E';
    end
    if center180
        dateChar = 'C';
    else
        dateChar = 'E';
    end
    gName = sprintf('D%s%04ixP%s%04i',dateChar,XYVec(1),poleChar,XYVec(2));
end
end