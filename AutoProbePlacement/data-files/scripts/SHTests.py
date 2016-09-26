 # -*- coding: utf-8 -*-
# from sympy.mpmath import *
from scipy import misc

from helper import *
import numpy
import math
import cPickle
import sys


height = 256;
width = 512;

maxBand = 9;
numCoeffs = 0;
for i in range(maxBand+1):
	numCoeffs = numCoeffs + (2*i + 1);

shArray = numpy.zeros((numCoeffs, height,width));
iArray = numpy.zeros((height,width));
jArray = numpy.zeros((height,width));
uArray = numpy.zeros((height,width));
vArray = numpy.zeros((height,width));
xArray = numpy.zeros((height,width));
yArray = numpy.zeros((height,width));
zArray = numpy.zeros((height,width));
longitudeArray = numpy.zeros((height,width));
latitudeArray = numpy.zeros((height,width));
epsilon = 0.01;
total = height * width * numCoeffs;

perten = math.floor(total / 10);
percent = math.floor(total / 100);
progress = 0;

for i in range(height):
	print(repr(i)+"/"+repr(height)+ " complete")
	for j in range(width):

		(u,v) = ijToUV_yup(i,j, height, width);

		(phi, theta) = equiUVtoPT_yup(u,v);
		(x,y,z) = sphericalToCartesian_yup(1,phi, theta);

		

		# if (i == 128):
		# 	print("i="+repr(i)+" j="+repr(j)+" u="+repr(u)+" v="+repr(v)+" phi="+repr(phi)+" theta="+repr(theta)+" x="+repr(x)+" y="+repr(y)+" z="+repr(z));

		# iArray[i][j] = float(i)/float(height);
		# jArray[i][j] = float(j)/float(width);
		# uArray[i][j] = (u+1)/2; #as u is between -1 and 1
		# vArray[i][j] = (v+1.0)/2.0; #as u is between -1 and 1
		# xArray[i][j] = x; #as u is between -1 and 1
		# yArray[i][j] = y; #as u is between -1 and 1
		# zArray[i][j] = z; #as u is between -1 and 1
		
		# if (phi == 0):
		# 	longitudeArray[i][j] = 1.0;
		# if (abs(phi-(pi/4)) < epsilon):
		# 	longitudeArray[i][j] = 1.0;
		# if (abs(phi-(pi/2)) < epsilon):
		# 	longitudeArray[i][j] = 1.0;
		# if (abs(phi-(3*pi/2)) < epsilon):
		# 	longitudeArray[i][j] = 1.0;
		# if (abs(phi+(pi/4)) < epsilon):
		# 	longitudeArray[i][j] = 1.0;
		# if (abs(phi+(pi/2)) < epsilon):
		# 	longitudeArray[i][j] = 1.0;
		# if (abs(phi+(3*pi/2)) < epsilon):
		# 	longitudeArray[i][j] = 1.0;
		#  #as u is between -1 and 1
		# if (abs(theta-(pi/4)) < epsilon):
		# 	latitudeArray[i][j] = 1.0;
		# if (abs(theta-(3*pi/4)) < epsilon):
		# 	latitudeArray[i][j] = 1.0;	
		# if (theta < epsilon):
		# 	latitudeArray[i][j] = 1.0;	
		for k in range(numCoeffs):
			shArray[k][i][j] = SHxyz_yup(kToLM(k), x,y,z);

for i in xrange(numCoeffs):
	cPickle.dump(shArray[i], open('../SH/512x256/values_'+repr(i)+'.p', 'wb')) 

# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/iArray.png",iArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/jArray.png",jArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/uArray.png",uArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/vArray.png",vArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/xArray.png",xArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/yArray.png",yArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/zArray.png",zArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/latitudeArray.png",latitudeArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/longitudeArray.png",longitudeArray);
