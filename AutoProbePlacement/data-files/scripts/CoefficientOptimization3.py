# -*- coding: utf-8 -*-
# Optimize coefficients for probes on a regular grid

import numpy as np;
from scipy.optimize import minimize;
from scipy.sparse import lil_matrix;
from scipy.sparse.linalg import inv
from scipy.sparse.linalg import lsmr;
from scipy.sparse.linalg import lsqr;
import helper
import platform;
import math;
np.set_printoptions(precision=5)
np.set_printoptions(suppress=True)
import gc
import time
import sys
import os
import random
LogPath = "../temp/PythonLog.txt";

global COLOR;
COLOR = -1;
#########################
### CONVENIENCE STUFF ###
#########################

def printSeparator():
	print("---------------------------");

# Take a list of strings and return their floating point representation
def LineToVectorF(line):
	vector = [float(x) for x in line.split()]
	# print(vector);
	return (vector);

# Take a list of ints and return their floating point representation
def LineToVectorI(line):
	vector = [int(x) for x in line.split()]
	# print(vector);
	return (vector);

# Append a vector to a list
def AppendToArray(array, vector):
	array.append(vector);
	
# Clear the ouput log.
def clearLog():
	writeStringToFile(LogPath, "", 'w');

# Print a Matrix with the prettiest formatting possible (pretty horrible imo)
def LogPrint(MatrixName, Matrix):
	try:
		writeStringToFile(LogPath, MatrixName+" \n" + assembleMatrixString(Matrix.tolist()));
	except:
		try:
			writeStringToFile(LogPath, MatrixName+" \n" + assembleMatrixString(Matrix));
		except:
			writeStringToFile(LogPath, MatrixName+" \n" + repr(Matrix) + '\n');

# Append a string to a file
def writeStringToFile(Path, String, option='a'):
	file = open(Path, option);
	file.write(String+"\n");
	file.close();
	
# Returns a string containing the pretty representation of a matrix
def assembleMatrixString(Matrix):

	String = "";

	for j in xrange(len(Matrix[0])):
		if (j == 0):
			String += "         "
		elif (j < 10):
			String += "          ";
		else:
			String += "         ";

		String += repr(j);
	String += "\n------";
	for j in xrange(len(Matrix[0])):
		String += "-----------";
	
	String += "\n";

	for i in xrange(len(Matrix)):
		if (i < 10):
			String += "0";
		String += repr(i) + " | ";
		for j in xrange(len(Matrix[0])):
			if ((Matrix[i][j] >= 0.0) and (repr(Matrix[i][j])[0] != "-")):
				String += " ";
			String += "%.3e" % Matrix[i][j] + " "
		String += "\n";
	
	return String;
	
######################
### ACTUAL PROGRAM ###
######################

# Load the information for each sample from the scene
def loadSamples():
	global ReferenceColors;
	samplesFile = open('../temp/samplePositions.txt','r');
	z = 1;

	for i in xrange(NumberOfReadSamples):
	
		# Load each line for this sample
		for j in xrange(NumberOfLinesPerSample):
			line = samplesFile.readline();

			while (line[0] == '#'):
				line = samplesFile.readline();

			if (j == 0): # Surface IDs
				# Put the relevant surface`s albedo in the list
				# SampleAlbedos.append(albedos[LineToVectorI(line)[0]]);
				pass
			if (j == 1): # Samples Positions
				SamplePositions.append(LineToVectorF(line));
			if (j == 3): # Normals
				SampleNormals.append(LineToVectorF(line));
			elif (j == 6): # Reference RGB
				ReferenceColors.append(LineToVectorF(line));

		samplesFile.readline(); #skip the empty line
	# Turn the ReferenceColors list-of-list in a matrix
	ReferenceColors = np.matrix(ReferenceColors);
	LogPrint("ReferenceColors", ReferenceColors);

def loadProbeCoefficients(probeStructure):
	global ProbeCoefficientsVector;
	for i in xrange(NumberOfProbes):
		print("\rPreparing probe " + repr(i) + " / " + repr(NumberOfProbes)),;

		for j in xrange(NumberOfCoefficientsPerProbe):
			ProbeCoefficientsVector[0].append(0);
			ProbeCoefficientsVector[1].append(0);
			ProbeCoefficientsVector[2].append(0);
	print("\n");

	for i in xrange(NumberOfProbes):
		print("\rLoading probe " + repr(i) + " / " + repr(NumberOfProbes)),;

		CoefficientFile = open(DataFilesPath+'ProbeStructures/'+probeStructure+'/coefficients/Coefficients_'+repr(i)+'.txt');
		# ProbeCoefficients[0].append([]);
		# ProbeCoefficients[1].append([]);
		# ProbeCoefficients[2].append([]);

		for j in xrange(NumberOfCoefficientsPerProbe):
			a = i*NumberOfCoefficientsPerProbe+j;
			ProbeCoefficientsVector[0][a] = float(CoefficientFile.readline());
			ProbeCoefficientsVector[1][a] = float(CoefficientFile.readline());
			ProbeCoefficientsVector[2][a] = float(CoefficientFile.readline());

	ProbeCoefficientsVector = np.array(ProbeCoefficientsVector);
	# Lolprint("ProbeCoefficientsVector.T", ProbeCoefficientsVector.T);
	print("\n");


# Save the 3 lists of coefficients by color to the coefficients file in "computed" ProbeStructure
def saveCoeffsToFile(R, G, B):
	global loadRandomProbes
	
	if (loadRandomProbes):
		structureID = "fromRandom_computed_"+repr(NumberOfReadSamples)+"samples_"+repr(NumberOfCoefficientsPerProbe)+"coeffs/";
	else:
		structureID = "computed_"+repr(NumberOfReadSamples)+"samples_"+repr(NumberOfCoefficientsPerProbe)+"coeffs/";

	savePath = "../ProbeStructures/"+structureID;
	if not (os.path.isdir(savePath)):
		 os.makedirs(savePath);
		 os.makedirs(savePath+"coefficients/");
		 
	for i in xrange(NumberOfProbes):
		file = open(savePath+"coefficients/Coefficients_"+repr(i)+".txt", 'w');
		for j in xrange(NumberOfCoefficientsPerProbe):
			file.write(repr(R[i*NumberOfCoefficientsPerProbe+j])+"\n");
			file.write(repr(G[i*NumberOfCoefficientsPerProbe+j])+"\n");
			file.write(repr(B[i*NumberOfCoefficientsPerProbe+j])+"\n");

	
# List of things we want to print to the log
def debugPrinting():
	global A;
	LogPrint("ReferenceColors",ReferenceColors[:,0]);
	LogPrint("np.dot(WeightMatrix, ProbeCoefficientsVector.T[:,0]).T", np.dot(WeightMatrix, ProbeCoefficientsVector.T[:,0]).T);
	LogPrint("np.dot(WeightMatrix, ProbeCoefficientsVector.T).T - ReferenceColors", np.dot(WeightMatrix, ProbeCoefficientsVector.T[:,0]).T - ReferenceColors[:,0]);
	LogPrint("np.sum(np.square(np.dot(WeightMatrix, ProbeCoefficientsVector.T[:,0]).T - ReferenceColors[:,0]))", np.sum(np.square(np.dot(WeightMatrix, ProbeCoefficientsVector.T[:,0]).T - ReferenceColors[:,0])));
			
# Obtain fitted coefficients by using WAY matrix`s pseudo inverse
def getFittedCoeffs():
	global A;
	
	# Gets fitted coefficients in one step
	A = np.dot(WeightMatrix.getI(), ReferenceColors).T;
	print(WeightMatrix.shape)
	print(ProbeCoefficientsVector.shape);
	print( ReferenceColors.shape)
	LogPrint("ReferenceColors.T[0].T", ReferenceColors.T[0].T);
	print( ReferenceColors.T[0].T.shape)

	# # print(np.dot(WeightMatrix, ProbeCoefficientsVector.T).shape);
	# A0 = lsqr(WeightMatrix, np.ndarray(ReferenceColors.T[0].T));
	# A1 = lsqr(WeightMatrix, np.ndarray(ReferenceColors.T[1].T));
	# A2 = lsqr(WeightMatrix, np.ndarray(ReferenceColors.T[2].T));

	# LogPrint("A", A);
	# Convert to list for easy outputting
	# A = A[0].T;
	# print A[0].shape
	
	B = A.tolist();
	
	# B0 = A0[0].tolist();
	# B1 = A1[0].tolist();
	# B2 = A2[0].tolist();
	# # LogPrint("A[0].shape", A.shape);
	# LogPrint("len(B)", len(B));
	# LogPrint("len(B[0])", len(B[0]));
	# LogPrint("Computed Coefficients from inversion (transposed)",A.T);
	# LogPrint("B", B);

	# LogPrint("B", B);
	# Output the optimized coefficients
	
	saveCoeffsToFile(B[0], B[1], B[2]);
	
	
	# for i in xrange(NumberOfProbes):
		# file = open(DataFilesPath+"ProbeStructures/computed/coefficients/Coefficients_"+repr(i)+".txt", 'w');
		# for j in xrange(NumberOfCoefficientsPerProbe):
			# file.write(repr(B[0][i*NumberOfCoefficientsPerProbe+j])+"\n");
			# file.write(repr(B[1][i*NumberOfCoefficientsPerProbe+j])+"\n");
			# file.write(repr(B[2][i*NumberOfCoefficientsPerProbe+j])+"\n");
			
# Get the lower corner for trilinear interpolation according to the position
def getSmallestCoordForInterpolation(position):
	x = math.floor((position[0] / step)) * step;
	y = math.floor((position[1] / step)) * step;
	z = math.floor((position[2] / step)) * step;
	
	return [x,y,z]

# Get 
def getInterpolationProbesCoords(position):
	smallestCoord = getSmallestCoordForInterpolation(position);
	currentCoord = smallestCoord;
	coords = [];
	coords.append(currentCoord);
	currentCoord = [smallestCoord[0] + step, smallestCoord[1], smallestCoord[2]];
	coords.append(currentCoord);
	currentCoord = [smallestCoord[0], smallestCoord[1] + step, smallestCoord[2]];
	coords.append(currentCoord);
	currentCoord = [smallestCoord[0], smallestCoord[1], smallestCoord[2] + step];
	coords.append(currentCoord);
	currentCoord = [smallestCoord[0] + step, smallestCoord[1], smallestCoord[2] + step];
	coords.append(currentCoord);
	currentCoord = [smallestCoord[0], smallestCoord[1] + step, smallestCoord[2] + step];
	coords.append(currentCoord);
	currentCoord = [smallestCoord[0] + step, smallestCoord[1] + step, smallestCoord[2]];
	coords.append(currentCoord);
	currentCoord = [smallestCoord[0] + step, smallestCoord[1] + step, smallestCoord[2] + step];
	coords.append(currentCoord);

	return coords;

# Assemble, line by line, the matrix holding all the info necessary for
# convolution of the samples
def fillWeightMatrix():
	global NumberOfReadSamples, NumberOfProbes, NumberOfCoefficientsPerProbe, WeightMatrix;
	global index_u, index_v;
	global SurfaceNormalsIJ

	for sampleIndex in xrange(NumberOfReadSamples):
		print("\rRead sample " + repr(sampleIndex+1) + " / " + repr(NumberOfReadSamples)),
		line = GlobalList[:];

		# Obtain the Sample Position
		position = SamplePositions[sampleIndex];
		
		# The amounts of probes in each direction
		repetitions = [21,11,21];

		# Get the positions of the probes to use for interpolation
		coords = getInterpolationProbesCoords(position);	

		# Get the importances of each maximal probes for interpolation (upper corners)
		ratios = [(position[0] - coords[0][0]) / step, (position[1] - coords[0][1]) / step, (position[2] - coords[0][2]) / step];

		# Specify the positions of the first probes in the array
		firstProbePosition = [-5,0,-5];
		lastProbePosition = [5,5,5];
		x = firstProbePosition[0];
		
		count = 0;
		weightSum = 0;
		
		# Iterate over all the interpolating probes
		for coord in coords:
			# print(("Interpolation probe coord: ", coord))
			badCoord = False;
			
			#Check if calculated interpolation position is outside level
			for i in xrange(3):
				if (coord[i] > lastProbePosition[i]):
					badCoord = True;
					break;
			# If it is move on to the next
			if (badCoord == True):
				count += 1;
				continue;
			
			# Set appropriate interpolation weight for this probe depending on its position in the list
			weight = 0;
			if (count == 0):
				weight = (1 - ratios[0])*(1-ratios[1])*(1-ratios[2]);
			elif (count == 1):
				weight = (ratios[0])*(1-ratios[1])*(1-ratios[2]); 
			elif (count == 2):
				weight = (1 - ratios[0])*(ratios[1])*(1-ratios[2]); 
			elif (count == 3):
				weight = (1 - ratios[0])*(1-ratios[1])*(ratios[2]); 
			elif (count == 4):
				weight = (ratios[0])*(1-ratios[1])*(ratios[2]); 
			elif (count == 5):
				weight = (1 - ratios[0])*(ratios[1])*(ratios[2]); 
			elif (count == 6):
				weight = (ratios[0])*(ratios[1])*(1-ratios[2]); 
			elif (count == 7):
				weight = (ratios[0])*(ratios[1])*(ratios[2]);

			weightSum += weight;

			# Get the index in the probes vector of the current probe
			index = int((2*coord[0] + 10) * repetitions[1] * repetitions[2] + \
					2 * coord[1] * repetitions[2] + \
					(2*coord[2] + 10)) * NumberOfCoefficientsPerProbe;

			# Iterate over all the SH bands
			lm = 0;
			for i in xrange(index, index + NumberOfCoefficientsPerProbe):
				(l,m) = helper.kToLM_new(lm);
				
				# Apply all multipliers for dot product
				factors = weight * helper.NDotOmegaCoeff(l) * helper.SHxyz_yup((l,m), SampleNormals[sampleIndex]);

				# print factors
				line[i] = factors;
				lm += 1;
							
			count += 1;

		WeightMatrix.append(line);
	
	# Get teh matrix from the list of lists
	# WeightMatrix = lil_matrix(WeightMatrix);
	WeightMatrix = np.matrix(WeightMatrix);

# Compute the cost function for one color using the relevant probe coefficient vector x
def computeRGBs(x):
	# How manyth iteration?
	global ReferenceColors
	global WeightMatrix
	global PassIndex
	global PRINTLOG
	
	# Global holding which color we`re looking at right now (kind of lazy but works)
	global COLOR

	# Dot product between WAY and probe coefficients - does all SH convolutions at once
	A = np.dot(WeightMatrix, x);
	
	# Get the colors by mulitplying element-to-element with surface albedos
	# Then compute the distances to the reference colors
	# B = A.T*SampleAlbedos[:][COLOR] - ReferenceColors[:,COLOR];
	B = A.T - ReferenceColors[:,COLOR];
	
	# Square all the results...
	maxdiff = 0;
	ind = 0
	C = np.square(B);


	# and sum them
	D = np.sum(C);

	# Control if we`re outputting to log
	if (PRINTLOG == 1):
		for i in xrange(len(C)):
			LogPrint("C["+repr(i)+"]",C[i])
			if (C[i] > maxdiff):
				maxdiff = C[i]
				ind = i
		LogPrint("np.dot(WeightMatrix, x).T", A.T);
		LogPrint("ReferenceColors[:,COLOR]", ReferenceColors[:,COLOR])
		LogPrint("squarediff", C);
	PassIndex += 1;
	
	# Return the cost!
	print(D);
	# exit()
	return D;

#################
# Program Start	#
#################

# Options
loadRandomProbes = False;
loadFittedProbes = False;
fitCoeffsByInversion = False;
fitCoeffsByMinimizer = False;
showCostFunctions = False;

if (sys.argv[1] == "fitI"):
	fitCoeffsByInversion = True;

elif (sys.argv[1] == "fitM"):
	loadFittedProbes = True;
	fitCoeffsByMinimizer = True;
	
elif (sys.argv[1] == "cost"):
	loadFittedProbes = True;
	showCostFunctions = True;

try:
	if (sys.argv[2] == "random"):
		loadRandomProbes = True;
		loadFittedProbes = False;

except:
	pass
# Clear log
clearLog();
if (platform.system() == 'Windows'):
	DataFilesPath = 'C:/libraries/g3d/samples/aTest/data-files/';
else:
	DataFilesPath = '/Users/joelpp/Documents/Maitrise/Code/G3D/G3D10/samples/aTest/data-files/';


# Initialize Parameters
NumberOfProbes = 4851;
NumberOfCoefficientsPerProbe = 4;
NumberOfElements = NumberOfProbes * NumberOfCoefficientsPerProbe;
NumberOfReadSamples = 1000;
NumberOfLinesPerSample = 7;
step = 0.5;

# Build the global list containing all 0's (for efficient creation of each line in the weight matrix)
GlobalList = [];
for i in xrange(NumberOfElements):
	GlobalList.append(0);

# Initialize lists we'll need
WeightMatrix = [];
ReferenceColors = [];
ProbeCoefficientsVector = [[],[],[]];
ProbeAssociations = [];
ProbeWeights = [];
SamplePositions = [];
SampleNormals = [];
SampleAlbedos = [];

PassIndex = 0;
PRINTLOG = 0;

print("Initializing...");

# Load the surface samples from file
loadSamples();
print("Finished loading Samples!");

# Load the probe Coefficients
if (loadFittedProbes):
	structureID = "computed_"+repr(NumberOfReadSamples)+"samples_"+repr(NumberOfCoefficientsPerProbe)+"coeffs/";
	loadProbeCoefficients(structureID);
elif (loadRandomProbes):
	loadProbeCoefficients("RandomCoefficients");
else:
	loadProbeCoefficients("RegularGrid");
	
print("Finished loading Probes!");

# Fill the weight matrix with appropriate interpolation weights / SH functions evaluations
fillWeightMatrix();
# WeightMatrix = [];
# for j in xrange(NumberOfReadSamples):
	# WeightMatrix.append([]);
	# for i in xrange(NumberOfElements):
		# WeightMatrix[j].append(random.uniform(0,1));
# WeightMatrix = np.matrix(WeightMatrix);
print("Filled WAY matrix 3!");

# Fit coefficients to reconstruct reference colors via matrix inversion
if (fitCoeffsByInversion):
	getFittedCoeffs();

# Compute cost functions with current probe coefficients
if (showCostFunctions):
	COLOR = 0;
	computeRGBs(ProbeCoefficientsVector.T[:,0])
	COLOR = 1
	computeRGBs(ProbeCoefficientsVector.T[:,1])
	COLOR = 2
	computeRGBs(ProbeCoefficientsVector.T[:,2])


if (fitCoeffsByMinimizer):
	COLOR = 0;
	resR = minimize(computeRGBs, ProbeCoefficientsVector.T[:,0], method='nelder-mead',options={'xtol': 1e-6, 'disp': True});
	COLOR = 1;
	resG = minimize(computeRGBs, ProbeCoefficientsVector.T[:,1], method='nelder-mead',options={'xtol': 1e-6, 'disp': True});
	COLOR = 2;
	resB = minimize(computeRGBs, ProbeCoefficientsVector.T[:,2], method='nelder-mead',options={'xtol': 1e-6, 'disp': True});
		
	LogPrint("resR",resR);
	LogPrint("resG",resG);
	LogPrint("resB",resB);

	saveCoeffsToFile(resR.x, resB.x, resG.x);
