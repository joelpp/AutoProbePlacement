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
import imageio
import helper;

# import matplotlib.pyplot as plt
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
	
def loadNormal(texel):
	# print(texel);
	normal = [x * 2 for x in texel];
	normal = [a - b for a, b in zip(normal, [1,1,1])];
	# print(normal);
	return normal;

def loadDepth(encodedDepth):
	return (1.0 - encodedDepth) * 10;
	
def loadPosition(texel):
	position = [x * 2 for x in texel];
	position = [a - b for a, b in zip(position, [1,1,1])];
	position = [x * 700 for x in position];

	return position;

def extractCoefficients(irradianceImage, normalImage, positionImage, depthImage, probePosition):
	height = len(irradianceImage);
	width = len(irradianceImage[0]);

	LSloan = [];
	LRam = [];
	coeffGrad = [];
	maxBand = 2;
	numCoeffs = 0;
	for i in range(maxBand+1):
		numCoeffs = numCoeffs + (2*i + 1);
	# numCoeffs = 2 * maxBand + 1;
	for i in range(numCoeffs):
		LSloan.append([0,0,0]);
		LRam.append([0,0,0]);
		coeffGrad.append([[0, 0, 0], [0, 0, 0], [0, 0, 0]]);

	pi = math.pi;
	total = height*width;
	perten = math.floor(total / 10);
	percent = math.floor(total / 100);
	progress = 0;
	domegasum = 0;
	domegaimg = numpy.zeros((height,width));
	lol = 0;
	logfile = open("C:/temp/log.txt", 'w');

	print("percent:"+repr(percent));
	# Loop over every pixel in the radiance map
	for h in xrange(height):
		for w in xrange(width):

			# Update progress bar if necessary
			if ((h*width+w) % perten == 0):
				print(repr(progress)+"% complete\r");
				progress += 10;

			
			(u,v) = ijToUV_yup(h,w, height, width);
			(phi, theta) = equiUVtoPT_yup(u,v);

			texelWSVector = sphericalToCartesian_yup(1,phi, theta);
			
			# Calculate differential solid angle element. Not sure this is right.
			domega = (pi/height)*(2*pi/width)*sin(theta);
			domegaimg[h][w] = domega;
			domegasum += domega;
			# logfile.write(repr((h,w,u,v,phi,theta,domega)) + "\n");

			radiance = [helper.sRGBToRGBVal(x /255.0) for x in irradianceImage[h][w]];
			# print(radiance)
			if ((radiance[0] == 0) and (radiance[1] == 0) and (radiance[2] == 0)):
				continue;
			# normal = loadNormal(normalImage[h][w]);
			# position = loadPosition(positionImage[h][w]);

			# lol += 1;
			# s = [a - b for a, b in zip(position, probePosition)];
			# sNorm = vectorNorm(s);
			
			# directionVector = normalizeVector(s);
			
			# normalDotDirectionVector = dot(normal, directionVector);
			
			# dA = domega * sNorm * sNorm / normalDotDirectionVector;
			# g = normalDotDirectionVector / (sNorm ** 2);
			
			# gGradTerm0 = [(x / (sNorm ** 3)) for x in normal];
			# gGradTerm1Factor = 3 * dot(normal, s) / (sNorm ** 5);
			# gGradTerm1 = [x * gGradTerm1Factor for x in s];
			
			# gGrad = [a - b for a, b in zip(gGradTerm0, gGradTerm1)];
			# (dYdx, dYdy, dYdz) = computeSHGradients(normal);

			# # dYdx = computedYdx(normal);
			# # dYdy = computedYdy(normal);
			# # dYdz = computedYdz(normal);
			# # # Step up in the integral for the 9 coefficients * 3 colors.
			for k in range(numCoeffs):
				# shFunctionEvaluation = SHxyz_yup(kToLM(k), tuple(directionVector));

				# shFunctionEvaluationGrad = [dYdx[k], dYdy[k], dYdz[k]];
				
				# integralFactor0Term0 = [x * g for x in shFunctionEvaluationGrad]
				# integralFactor0Term1 = [x * shFunctionEvaluation for x in gGrad]
				
				# integralFactor0 = [a + b for a, b in zip(integralFactor0Term0, integralFactor0Term1)]
				
				# integralTerm = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
				# integralTerm[0] = [x * radiance[0] * dA for x in integralFactor0];
				# integralTerm[1] = [x * radiance[1] * dA for x in integralFactor0];
				# integralTerm[2] = [x * radiance[2] * dA for x in integralFactor0];
				
				# for channel in xrange(3):
				# 	coeffGrad[k][channel] = [(a + b) * 1 for a,b in zip(coeffGrad[k][channel], integralTerm[channel])]
				
				# if ((k == 0) and (h == 15) and (w == 0)):
				# 	ss = "";
				# 	ss += "DEBUG INFO ";
				# 	ss += repr(position) + "\n";
				# 	ss += repr(normal) + "\n";
				# 	ss += repr(s) + "\n";
				# 	ss += repr(normalDotDirectionVector) + "\n";
				# 	ss += repr(gGrad) + "\n";
				# 	ss += repr(dYdx[0]) + "\n";
				# 	ss += repr(shFunctionEvaluationGrad) + "\n";
				# 	ss += repr(integralFactor0) + "\n";
				# 	ss += repr(integralTerm[0]) + "\n";
				# 	ss += "END DEBUG INFO \n";
				# 	print(ss)
				
				shFunctionEvaluationForTexel = SHxyz_yup(kToLM(k), texelWSVector);
				(l,m) = helper.kToLM(k);
				# if ((h ==5 ) and (w == 5)):
				# 	print((radiance, texelWSVector, domega, k, shFunctionEvaluationForTexel));
				LSloan[k][0] += domega * radiance[0] * helper.NDotOmegaCoeff(l) * shFunctionEvaluationForTexel;
				LSloan[k][1] += domega * radiance[1] * helper.NDotOmegaCoeff(l) * shFunctionEvaluationForTexel;
				LSloan[k][2] += domega * radiance[2] * helper.NDotOmegaCoeff(l) * shFunctionEvaluationForTexel;
				# if (k == 0):
					# logfile.write(repr(i) + " " + repr(j) + " " + repr(domega * radiance[0] * shFunctionEvaluationForTexel) + " " + repr(domega * radiance[1] * shFunctionEvaluationForTexel) + " " + repr(domega * radiance[2] * shFunctionEvaluationForTexel) + "\n");
					# logfile.write(repr(h) + " " + repr(w) + " " + repr(domega) + " " + repr(radiance) + " " + repr(shFunctionEvaluationForTexel) + " " + repr(domega * radiance[0] * shFunctionEvaluationForTexel) + " " + repr(domega * radiance[1] * shFunctionEvaluationForTexel) + " " + repr(domega * radiance[2] * shFunctionEvaluationForTexel) + "\n");
	print(domegasum);
	return (LSloan, coeffGrad);

#------------------------------
#|		PROGRAM START 		  |	  	
#------------------------------
sys.setrecursionlimit(10000);

offset = (0,0,0);
# expNumber = int(sys.argv[1]);

(sceneName, structureName) = helper.getSceneAndStructureFromCommandLineArguments(sys.argv);
globalInfo = helper.readSceneAndProbeStructureInfo(sceneName, structureName);
rootPath = "../Scenes/" + sceneName + "/ProbeStructures/" + structureName;

if (len(sys.argv) == 4):
	id = int(sys.argv[3]);
else:
	id = -1;
	
# probeList = open(rootPath+"/probeList.txt", "r");
probeCount = 0;
with open(rootPath+"/probeList.txt") as f:
	lines = f.read().splitlines()

for line in lines:

	if ((id != -1) and (id != probeCount)):
		probeCount += 1;
		continue;

	start = time.time();

	splitLine = line.split();
	probePosition = [float(x) / globalInfo["scale"] for x in splitLine];

	path = "/probes/Probe_"+repr(probeCount)+".png";
	# path = "/probes/Probe_"+repr(probeCount)+".exr";
	normalPath = "/Normals/Normal_"+repr(probeCount)+".exr";
	depthPath = "/Depths/Depth_"+repr(probeCount)+".exr";
	positionPath = "/Positions/Position_"+repr(probeCount)+".exr";

	irradianceImage = imageio.imread(rootPath + path)
	# normalImage = imageio.imread(rootPath + normalPath)
	# depthImage = 0;#imageio.imread(rootPath + depthPath)
	# #positionImage = misc.imread(rootPath + positionPath)
	# positionImage = imageio.imread(rootPath + positionPath)

	#------------------------------
	#|	QUADRATURES INTEGRATION  |	  	
	#------------------------------i=
	savePath = rootPath+"/coefficients/Coefficients_"+repr(probeCount);
	gradSavePath = rootPath+"/coefficients/CoefficientsGradients_"+repr(probeCount);
	# savePath = "/Users/joelpp/Documents/Maitrise/Images/grace/coefficients";
	# if (os.path.isfile(savePath+'.p')):
		# print("Probe #"+repr(probeCount)+" already Mapped! continuing...");
		# probeCount = probeCount + 1;
		# continue;
	
	print("Beginning Quadratures Integration for probe "+repr(probeCount));

	# (LSloan, coeffGrad) = extractCoefficients(irradianceImage, normalImage, positionImage, depthImage, probePosition);
	(LSloan, coeffGrad) = extractCoefficients(irradianceImage, 0, 0, 0, probePosition);

	#Dump resulting coefficients in their own file and the console
	cPickle.dump(LSloan, open(savePath+'.p', 'wb')) 

	output = open(savePath+'.txt', 'wb');
	for i in range(len(LSloan)):
		for j in range(len(LSloan[i])):
			output.write(repr(LSloan[i][j])+"\n");
	output.close();
	
	output = open(gradSavePath+'.txt', 'wb');
	for coeff in xrange(len(coeffGrad)):
		output.write("// Coeff #" + repr(coeff) + "\n" );
		for channel in xrange(len(coeffGrad[coeff])):
			output.write("// Channel " + intToChannelString(channel) + "\n" );

			for axis in xrange(len(coeffGrad[coeff][channel])):
				output.write( repr( coeffGrad[coeff][channel][axis] * 1 ) + "\n" );
	output.close();
	
	
	# exit();
	print("L final vector for Quadratures (Sloan): "+repr(LSloan));
	probeCount = probeCount + 1;
	stop = time.time();
	print("time taken: "+repr(stop-start)+"s");
		# Visual represntation of the dw's.
		# misc.imsave("/Users/joelpp/Documents/Maitrise/Images/IrradianceMaps/domega.png", domegaimg);
		# show("domegasum",domegasum);
	# exit();


