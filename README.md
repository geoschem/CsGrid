CSGrid repository README

Seb Eastham created this MATLAB library of cubed sphere grid functions. 
Below is a Crash-course in usage provided by Seb. Additional things to
note include:
     (1) Not all tilefiles are included in the repository due to file size.
         Only the cubed sphere <-> 2x25 tile file is currently included.
     (2) Some example scripts are provided in the exampleScripts/
         subdirectory. 
     (3) This repository is still under development.
     (4) Contact Seb Eastham or the GEOS-Chem Support Team for additional
         files or for more information.
	 
~Lizzie Lundgren, 9/13/2016
_____________________________________________________________________________

1) EXTRACT VARIABLE FROM GCHP OUTPUT

Extract a variable from a GCHP output file, say 'testData'. For the moment, 
generate an arbitrary one with the dimensions of a C24 simulation, e.g.

 >    testData = zeros(24,144,72);
 >    for iFace=1:6
 >    for iLayer = 1:72
 >    testData(:,(1:24)+(iFace-1)*24,iLayer) = iFace*iLayer;
 >    end
 >    end

Note that the dimensions of a cubed sphere file with side length NX, NZ layers,
and T samples are [NX,NX*6,NZ,T].

_____________________________________________________________________________

2) PLOT VARIABLE IN NATIVE FORMAT (C24)

Plot the data in its native format using the plotCSLayer function. For example:

 > plotCSLayer(testData(:,:,1),'projection','flat');

_____________________________________________________________________________

3) REGRID VARIABLE FROM C24 TO 2X25

There's no sensible way to plot the global variation of quantities with 
height in a cubed-sphere grid other than to do some regridding. To regrid from 
C24 to 2x2.5, run:

 > LLData = regridConservative(testData,genGridSpec('gmao2x25','geos5'));

_____________________________________________________________________________

4) PLOT THE REGRIDDED VARIABLE

You can then plot it either using your preferred functions or using those that 
I've packaged here, e.g.

 > plotGrid(LLData,'zonal');

_____________________________________________________________________________

5) EXTRACT ONE LAYER OF THE PLOT

To get one layer of the regridded data, run:

 > plotGrid(LLData(:,:,1),'layer');

_____________________________________________________________________________

6) COMPARE GRID AREAS (in the lat-lon and CS frameworks) 

To get lat-lon statistics, run:

 > gSpec = genGridSpec('gmao2x25','geos5'); 

This will generate a file with a custom "grid specification" class (used 
implicitly in the regridding, earlier). You can then get the latitude edges, 
longitude edges and grid cell areas from:

 > gSpec.latEdge
 > gSpec.lonEdge
 > gSpec.gridArea

To get the equivalent data for the cubed-sphere, run:

 > [lonEdge,latEdge] = calcCSGrid(24); 
 > gridArea = calcCSArea(lonEdge,latEdge); 

The total area from the CS grid and the lat-lon grid should be equal to 
within some very large number of decimal places.

_____________________________________________________________________________

There's an awful lot of extra functionality which I haven't really gone 
into here in detail (e.g. you can update a CS layer plot much more quickly 
after the first plot, using the function "updateCSPatch" with the output 
from plotCSLayer; this can save a lot of time when you want to see how a 
variable evolves or even if you just want to look at lots of variables 
one after the other).

Regards,

Seb