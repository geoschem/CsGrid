function [ outData ] = processGCHPData( inData )
%PROCESSGCHPDATA Converts GCHP data to something more user-friendly
%   Incoming data is assumed to be in the following format:
%       Ix6IxLxT
%   with the vertical co-ordinate, L, reversed. Outgoing data will be:
%       IxIx6xLxT

inSize = size(inData);
assert(inSize(2) == 6*inSize(1),'processGCHPData:badSize','Data is not correctly sized');
nPerSide = inSize(1);
if length(inSize) >= 3
    nLev = inSize(3);
    if length(inSize) >= 4
        nSamples = inSize(4);
        assert(length(inSize) == 4,'processGCHPData:extraDims','Too many dimensions');
    else
        nSamples = 1;
    end
else
    nLev = 1;
end

outSize = [nPerSide,nPerSide,6,nLev,nSamples];
outData = zeros(outSize,'like',inData);

for iSample = 1:nSamples
    for iLev = 1:nLev
        inLev = nLev + 1 - iLev;
        outData(:,:,:,iLev,iSample) = reorderCS(inData(:,:,inLev,iSample),true);
    end
end

end

