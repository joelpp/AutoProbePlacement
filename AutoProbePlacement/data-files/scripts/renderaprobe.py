from scipy import misc
import numpy
import cPickle
import math
import time;
import os;
import sys
import helper
#------------------------------
#|		USEFUL FUNCTIONS	  |	  	
#------------------------------


#------------------------------
#|		PROGRAM START 		  |	  	
#------------------------------
start = time.time();
(sceneName, structureName) = helper.getSceneAndStructureFromCommandLineArguments(sys.argv);
globalInfo = helper.readSceneAndProbeStructureInfo(sceneName, structureName);
rootPath = "../Scenes/" + sceneName + "/ProbeStructures/" + structureName;
probeNumber = 0;

with open(rootPath+"/probeList.txt") as f:
	lines = f.read().splitlines()

for line in lines:

	height = 256;
	width = 512;

	maxBand = 2;
	numCoeffs = 0;
	nullList = [];
	numChannels = 3;

	# for i in range(maxBand+1):
		# numCoeffs = numCoeffs + (2*i + 1);
	numCoeffs = 9
	for i in xrange(numCoeffs):
		nullList.append([0,0,0]);
	phongExponent = 1;

	offset = (0,0,0);

	savePath = rootPath + "/renders/Render_"+repr(probeNumber)+".png";
		

	L = cPickle.load(open(rootPath + '/coefficients/Coefficients_'+repr(probeNumber)+'.p', 'rb'));

	print("Probe #"+repr(probeNumber));
	print(type(L));
	print(L);

		#Initialize image arrays.
	RGB = numpy.zeros((height, width, numChannels));


	if (L == nullList):
		print("Coeffs are empty. Continuing...");
		misc.imsave(savePath, RGB);
		continue;

		#Initialize some values for the counter.
	total = height*width;
	perten = math.floor(total / 10);
	progress = 0;
	epsilon = 0.0001;

		# Loop over all them texels.
	for i in range(height):
		for j in range(width):

			#Increment progress counter if appropriate
			if ((i*width+j) % perten == 0):
				print(repr(progress)+"% complete")
				progress += 10;

			#Convert from i,j to u,v coordinates
			# (u,v) = ijToUV(i,j, height, width);
			(u,v) = helper.ijToUV_yup(i,j, height, width);

			#Convert from u,v to theta,phi coordinates
			# (phi, theta) = equiUVtoPT(u,v);
			(phi, theta) = helper.equiUVtoPT_yup(u,v);

			#Convert from theta,phi to x,y,z coordinates
			# (x,y,z) = sphericalToCartesian(1,phi, theta);
			(x,y,z) = helper.sphericalToCartesian_yup(1,phi, theta);

		
			# Color the point according to the SH coefficients as described in Ramamoorthi's article.

			if (numChannels == 4):
				RGB[i][j][3] = 1.0;
			for k in range(3):
				# This is the rendering equation using Sloan's coefficients
				# RGB[i][j][k] = 0.4288 * L[8][k] * (x*x - y*y) + 0.247 * 3 * L[6][k] * z*z + 2.7834 * L[0][k] - 0.247 * L[6][k] + 2.6943 * ( L[4][k]*x*y - L[7][k]*x*z - L[5][k]*y*z) + 3.2149 * (-L[3][k] * x - L[1][k]*y + L[2][k]*z);
				
				for coeff in range(numCoeffs):
		
					# RGB[i][j][k] += L[coeff][k] * phongCoeffs(kToLM(coeff), 1) * SHxyz(kToLM(coeff), x,y,z);
					# RGB[i][j][k] += L[coeff][k] * phongCoeffs(kToLM(coeff), phongExponent) * SHBetterLookup(i,j,coeff);
					# print(L[coeff][k])
					# print(helper.NDotOmegaCoeff(coeff))
					# print(helper.SHxyz_yup(helper.kToLM(coeff),(x,y,z)))
					(l,m) = helper.kToLM(coeff);
					RGB[i][j][k] += L[coeff][k] * helper.NDotOmegaCoeff(l) * helper.SHxyz_yup((l,m),(x,y,z));


		

	misc.imsave(savePath, RGB);
	probeNumber += 1;

# knock knock. Time to go!
end = time.time();
print("Elapsed: "+repr((end-start)));
