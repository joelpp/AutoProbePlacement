from sympy.mpmath import *
from scipy import misc

from helper import *
import numpy
import cPickle
import math

path = "Mitsuba/RedGreen/SphericalCamera_BetterScaleObj_0_0_0.png";
img = misc.imread("/Users/joelpp/Documents/Maitrise/Images/"+path);

# height = shArray.shape[0];
# width = shArray.shape[1];

height = 128;
width = 256;



maxBand = 9;
numCoeffs = 0;
for i in range(maxBand+1):
	numCoeffs = numCoeffs + (2*i + 1);
total = height*width * numCoeffs;


perten = math.floor(total / 10);
percent = math.floor(total / 100);
progress = 0;
quarter = width / 4;
arr = numpy.zeros((height,width));

for i in xrange(numCoeffs):
	for j in range(height):
		for k in range(width):
			arr[j][k] = (SHBetterLookup_small(j,k,i)+1)/2.0;
			# arr[j][(k)%width] = (shArray[j][k][i]+1)/2.0;
			# if (i ==0 ):
				# arr[j][(k)%width] = (0.2820947+1)/2.0;
	misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/256x128/SHRender"+repr(i)+".png", arr);
	


# for i in range(height):
# 	for j in range(width):

# 		# Update progress bar if necessary
# 		if ((i*width+j) % percent == 0):
# 			print(repr(progress)+"% complete\r");
# 			progress += 1;

# 		# print((i*width+j) % percent);
# 		#Convert from i,j to u,v coordinates
# 		(u,v) = ijToUV(i,j, height, width);
# 		# print("i: "+repr(i)+"j: "+repr(j)+"u: "+repr(u)+"v: "+repr(v));
# 		#Convert from u,v to theta,phi coordinates
# 		(phi, theta) = equiUVtoPT(u,v);

# 		#This is taken from prefilter.c and temporary
# 		# R = math.sqrt(u*u+v*v);  
# 		# phi = math.pi*R;
# 		# theta = math.atan2(u,v);

# 		#Convert from theta,phi to x,y,z coordinates
# 		(x,y,z) = sphericalToCartesian(1,phi, theta);

# 		for k in range(numCoeffs):
# 			lmtuple = kToLMLegendre(k);
# 			legendreArray[i][j][k] = legenp(lmtuple[0], lmtuple[1], math.cos(theta));

# cPickle.dump(legendreArray, open('/Users/joelpp/Documents/Maitrise/Code/IrradianceMapping/LegendrePolynomials/polynomials.p', 'wb')) 
