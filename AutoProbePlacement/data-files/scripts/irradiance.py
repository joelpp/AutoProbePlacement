from __future__ import division
from __future__ import print_function
import math
import time
import random
import numpy
import cPickle
import sys
import struct
from scipy import misc
from scipy import ndimage
from helper import *
import os.path
# import matplotlib.pyplot as plt



def extractCoefficients(img):
	height = len(img);
	width = len(img[0]);

	LSloan = [];
	LRam = [];

	maxBand = 2;
	numCoeffs = 0;
	for i in range(maxBand+1):
		numCoeffs = numCoeffs + (2*i + 1);
	# numCoeffs = 2 * maxBand + 1;
	for i in range(numCoeffs):
		LSloan.append([0,0,0]);
		LRam.append([0,0,0]);

	pi = math.pi;
	total = height*width;
	perten = math.floor(total / 10);
	percent = math.floor(total / 100);
	progress = 0;
	domegasum = 0;
	domegaimg = numpy.zeros((height,width));
	print("percent:"+repr(percent));
	print(numCoeffs);
	# Loop over every pixel in the radiance map
	for i in range(height):
		for j in range(width):

			# Update progress bar if necessary
			if ((i*width+j) % perten == 0):
				print(repr(progress)+"% complete\r");
				progress += 10;

			# print((i*width+j) % percent);
			#Convert from i,j to u,v coordinates
			#
			# (u,v) = ijToUV(i,j, height, width);
			(u,v) = ijToUV_yup(i,j, height, width);

			# print("i: "+repr(i)+"j: "+repr(j)+"u: "+repr(u)+"v: "+repr(v));
			#Convert from u,v to theta,phi coordinates
			#
			# (phi, theta) = equiUVtoPT(u,v);
			(phi, theta) = equiUVtoPT_yup(u,v);

			#This is taken from prefilter.c and temporary
			# R = math.sqrt(u*u+v*v);  
			# phi = math.pi*R;
			# theta = math.atan2(u,v);

			#Convert from theta,phi to x,y,z coordinates
			#
			# (x,y,z) = sphericalToCartesian(1,phi, theta);
			texelWSVector = sphericalToCartesian_yup(1,phi, theta);
			
			# Calculate differential solid angle element. Not sure this is right.
			domega = (pi/height)*(2*pi/width)*sin(theta);
			domegaimg[i][j] = domega;
			domegasum += domega;

			# print("SH(1,1): "+repr(SH((1,0), x, y, z)));
			# print("SHtp(1,1): "+repr(SHtp((1,0), phi, theta)));

			# For use with GRACE_PROBE ONLY
			if ((img[i][j][0] == 0) and (img[i][j][1] == 0) and (img[i][j][2] == 0)):
				continue;

			# r = struct.unpack('f', f.read(4));
			# g = struct.unpack('f', f.read(4));
			# b = struct.unpack('f', f.read(4));

			# a[i][j][0] = random.random();
			# a[i][j][1] = random.random();
			# a[i][j][2] = random.random();
			# if ((r == 0) and (g == 0) and (b == 0)):
			# 	continue;

			# # # Step up in the integral for the 9 coefficients * 3 colors.
			for k in range(numCoeffs):
				shFunctionEvaluationForTexel = SHxyz_yup(kToLM(k), texelWSVector);
				
				LSloan[k][0] += domega * float(img[i][j][0])/255 * shFunctionEvaluationForTexel;
				LSloan[k][1] += domega * float(img[i][j][1])/255 * shFunctionEvaluationForTexel;
				LSloan[k][2] += domega * float(img[i][j][2])/255 * shFunctionEvaluationForTexel;

				# LSloan[k][0] += domega * float(img[i][j][0])/255 * SHBetterLookup(i,j,k);
				# LSloan[k][1] += domega * float(img[i][j][1])/255 * SHBetterLookup(i,j,k);
				# LSloan[k][2] += domega * float(img[i][j][2])/255 * SHBetterLookup(i,j,k);
				# LSloan[k][3] = 0;
	print(LSloan)
	exit()
	return LSloan;

#------------------------------
#|		PROGRAM START 		  |	  	
#------------------------------
sys.setrecursionlimit(10000);

offset = (0,0,0);
# expNumber = int(sys.argv[1]);

sceneName = "2Lights"
structureName = "1ProbeAtOrigin"
rootPath = "../Scenes/" + sceneName + "/ProbeStructures/" + structureName;


probeList = open(rootPath+"/probeList.txt", "r");
probeCount = 0;
for line in probeList:
	# if (probeCount < 150):
		# probeCount += 1;
		# continue
# for i in xrange(1):
	start = time.time();

	# splitLine = line.split();

	path = "/probes/Probe_"+repr(probeCount)+".png";

	img = misc.imread(rootPath+path)
	# img = misc.imread("/Users/joelpp/Documents/Maitrise/Images/grace/grace_probe.png");
	#------------------------------
	#|	QUADRATURES INTEGRATION  |	  	
	#------------------------------i=
	savePath = rootPath+"/coefficients/Coefficients_"+repr(probeCount);
	# savePath = "/Users/joelpp/Documents/Maitrise/Images/grace/coefficients";
	# if (os.path.isfile(savePath+'.p')):
		# print("Probe #"+repr(probeCount)+" already Mapped! continuing...");
		# probeCount = probeCount + 1;
		# continue;
	
	print("Beginning Quadratures Integration for probe "+repr(probeCount));

	LSloan = extractCoefficients(img);

	#Dump resulting coefficients in their own file and the console
	cPickle.dump(LSloan, open(savePath+'.p', 'wb')) 
	# cPickle.dump(LRam, open('QuadratureCoefficientsRamamoorthi_'+repr(s)+'_'+repr(t)+'.p', 'wb')) 

	output = open(savePath+'.txt', 'wb');
	for i in range(len(LSloan)):
		for j in range(len(LSloan[i])):
			output.write(repr(LSloan[i][j])+"\n");

	print("L final vector for Quadratures (Sloan): "+repr(LSloan));
	probeCount = probeCount + 1;
	stop = time.time();
	print("time taken: "+repr(stop-start)+"s");
		# Visual represntation of the dw's.
		# misc.imsave("/Users/joelpp/Documents/Maitrise/Images/IrradianceMaps/domega.png", domegaimg);
		# show("domegasum",domegasum);



