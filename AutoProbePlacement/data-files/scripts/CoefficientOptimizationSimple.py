# -*- coding: utf-8 -*-
# Optimize coefficients for probes on a regular grid

import numpy as np;
from scipy.optimize import minimize, fmin, fmin_bfgs;
from scipy.sparse import lil_matrix, bsr_matrix, csr_matrix;
from scipy.sparse.linalg import inv, spsolve;
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
RED = 0;
GREEN = 1;
BLUE = 2;


DFrequencies = [];
Colors = [];
for i in xrange(3000):
	DFrequencies.append(0);

class Scene:
	def __init__(self, name, firstProbePosition, lastProbePosition, dimensions, step, scale):
		self.name = name;
		self.dimensions = dimensions;
		self.firstProbePosition = firstProbePosition;
		self.lastProbePosition = lastProbePosition;
		self.step = step;
		self.scale = scale;

def randomWeights():
	return [[0.8343103946030044, 0.685366951997276, 0.1412105940974021, 0.640776821605616, 0.5297586346926327, 0.4865722676612464, 0.4657218157255085, 0.8409810907359956, 0.7707315803804138], \
			[0.47651812547795713, 0.5973510865341789, 0.7339525638567312, 0.6306252914171964, 0.7398008814658011, 0.579932420848727, 0.1853927638444296, 0.3875308412073958, 0.040692633755435526], \
			[0.1727596643864473, 0.48651451785483124, 0.9030659976022873, 0.08216917814707903, 0.27258519316516916, 0.8664425377147795, 0.8617384641123248, 0.7124378475328639, 0.7521652778748018], \
			[0.059853148310527904, 0.2901058652845424, 0.5017090699975343, 0.9635157544378622, 0.24317488706534662, 0.40659184111519053, 0.20073498570283865, 0.29904574632004477, 0.2956003374285453], \
			[0.5032137397383104, 0.11597179462815221, 0.5032798715192168, 0.06895927670139712, 0.5877727114657686, 0.46134689269125706, 0.14293084713437842, 0.8937383095521069, 0.3761742868576353], \
			[0.49609460842443076, 0.4479032703830026, 0.5164817361989202, 0.3220755663169278, 0.7211059850934786, 0.310201687303577, 0.7112224244435797, 0.9306228714651787, 0.6708467115946546], \
			[0.73618527535356, 0.350753731895979, 0.23665868042470373, 0.7411827760474671, 0.3153475235115749, 0.5249109105434353, 0.024763749855696293, 0.2523419333364374, 0.21060989340458658], \
			[0.07079092185489266, 0.8277069436172073, 0.5611716633679361, 0.9274978338046512, 0.08788158178183325, 0.15351759684942723, 0.08646905317175446, 0.3280683246521916, 0.25135815188583843], \
			[0.3977192292650714, 0.6522208543520609, 0.15150833580532197, 0.29692844399078466, 0.8553588738025545, 0.9532709079062245, 0.3235657760295175, 0.45618799629683693, 0.8150213045998811]];

def printSeparator():
	print("---------------------------");

# Take a list of strings and return their floating point representation
def LineToVectorF(line):
	vector = [float(x) for x in line.split()]
	# print(vector);
	return (vector);
	
def LineToVectorF(line, multiplier):
	vector = [(float(x) * multiplier) for x in line.split()]
	return (vector);
	
# Take a list of ints and return their floating point representation
def LineToVectorI(line):
	vector = [int(x) for x in line.split()]
	# print(vector);
	return (vector);

def TupleLineToVectorF(line):
	fixed = line.replace(',', '');
	fixed = fixed.replace('(', '');
	fixed = fixed.replace(')', '');
	fixed = fixed.replace('[', '');
	fixed = fixed.replace(']', '');
	fixed = fixed.split();
	fixed.pop(0);
	vector = [int(x) for x in fixed];
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

def LogPrintList(ListName, List):
	String = "";
	String += ListName;
	String += "\n------\n";
	
	for i in xrange(len(List)):
		String += repr(i) + " : " + repr(List[i]) + "\n"
	
	writeStringToFile(LogPath, String);
	
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
		elif (j < 100):
			String += "         ";		
		elif (j < 1000):
			String += "        ";
		elif (j < 10000):
			String += "       ";		
		elif (j < 100000):
			String += "      ";
			
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
def loadSamples(sampleSetName):
	global ReferenceColors, scene;
	samplesFile = open('../Scenes/'+scene.name+'/SampleSets/' + sampleSetName + '/SamplePositions.txt','r');
	irradianceFile = open('../Scenes/'+scene.name+'/SampleSets/' + sampleSetName + '/IrradianceResults2.txt', 'r');
	z = 1;

	# for i in xrange(NumberOfReadSamples * 1000):
	for i in xrange(NumberOfReadSamples):
	
		# if (i < 30000):
			# irradianceFile.readline();
			# for k in xrange(NumberOfLinesPerSample):
				# samplesFile.readline();
			# samplesFile.readline(); #skip the empty line
			# continue;
		# Initialize the sample info
		SamplePosition = [];
		SampleNormal = [];
		ReferenceColor = [];
		
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
				# SamplePositions.append(LineToVectorF(line, 1));
				SamplePosition = [i * scene.scale for i in LineToVectorF(line, 1)]; 
			if (j == 4): # Normals
				# SampleNormals.append(LineToVectorF(line, 1));
				SampleNormal = LineToVectorF(line, 1);
			elif (j == 6): # Reference RGB
				# ReferenceColors.append(LineToVectorF(line));
				irrLine = irradianceFile.readline();
				# ReferenceColors.append(LineToVectorF(line));
				# ReferenceColors.append(LineToVectorF(irrLine, 1. / 255.));
				ReferenceColor = LineToVectorF(irrLine, 1.);

		samplesFile.readline(); #skip the empty line
		
		# if (ReferenceColor == [0,0,0]):
			# continue;
		# if ((SamplePosition[1] == 0) and (SampleNormal[1] == -1)):
				# ReferenceColor = [0,0,0]
		
		# Append the gathered info to the relevant buffers
		SamplePositions.append(SamplePosition);
		SampleNormals.append(SampleNormal);
		ReferenceColors.append(ReferenceColor);
		
		if (len(SamplePositions) == NumberOfReadSamples):
			break;
		
	# Turn the ReferenceColors list-of-list in a matrix
	# ReferenceColors = np.matrix(ReferenceColors);
	# LogPrint("ReferenceColors", ReferenceColors);
	# LogPrint("SamplePositions", SamplePositions);
	# LogPrint("SampleNormals", SampleNormals);

def loadProbeCoefficients(probeStructure):
	global ProbeCoefficientsVector;
	for i in xrange(NumberOfProbes):
		for j in xrange(NumberOfCoefficientsPerProbe):
			ProbeCoefficientsVector[0].append(0);
			ProbeCoefficientsVector[1].append(0);
			ProbeCoefficientsVector[2].append(0);
	print("\n");

	for i in xrange(NumberOfProbes):
		print("\rLoading probe " + repr(i) + " / " + repr(NumberOfProbes)),;

		CoefficientFile = open("../Scenes/" + scene.name + '/ProbeStructures/' + probeStructure + '/coefficients/Coefficients_' + repr(i) + '.txt');
		
		for j in xrange(NumberOfCoefficientsPerProbe):
			a = i*NumberOfCoefficientsPerProbe+j;
			ProbeCoefficientsVector[RED][a] = float(CoefficientFile.readline());
			ProbeCoefficientsVector[GREEN][a] = float(CoefficientFile.readline());
			ProbeCoefficientsVector[BLUE][a] = float(CoefficientFile.readline());

	ProbeCoefficientsVector = np.array(ProbeCoefficientsVector);
	# Lolprint("ProbeCoefficientsVector.T", ProbeCoefficientsVector.T);
	print("\n");


# Save the 3 lists of coefficients by color to the coefficients file in "computed" ProbeStructure
def saveCoeffsToFile(R, G, B):
	global loadRandomProbes, NumberOfReadSamples;
	global scene;

	StructureName = helper.buildStructurePath("Fitted_" + repr(NumberOfReadSamples) + "_Samples", scene.step, scene.dimensions, scene.firstProbePosition);
	StructurePath = "../Scenes/" + scene.name + "/ProbeStructures/" + StructureName;
	
	if not (os.path.isdir(StructurePath)):
		 os.makedirs(StructurePath);
		 os.makedirs(StructurePath+"/coefficients/");
		 
	infoFile = open(StructurePath + "/info.txt", "w");
	infoFile.write("step " + repr(scene.step) + "\n");
	infoFile.write("dimensions " + repr(scene.dimensions[0]) + " " +  repr(scene.dimensions[1]) + " " + repr(scene.dimensions[2]) + "\n");
	infoFile.write("firstProbePosition " + repr(scene.firstProbePosition[0]) + " " +  repr(scene.firstProbePosition[1]) + " " + repr(scene.firstProbePosition[2]) + "\n");

	probeListFile = open(StructurePath + "/probeList.txt", "w");
	
	
	for i in xrange(NumberOfProbes):
		
		probePosition = helper.positionFromIndex(i, scene.dimensions, scene.firstProbePosition);
		probeListFile.write(repr(probePosition[0]) + " " + repr(probePosition[1]) + " " + repr(probePosition[2]) + "\n");
	
		file = open(StructurePath+"/coefficients/Coefficients_"+repr(i)+".txt", 'w');
		for j in xrange(NumberOfCoefficientsPerProbe):
			file.write(repr(R[i*NumberOfCoefficientsPerProbe+j])+"\n");
			file.write(repr(G[i*NumberOfCoefficientsPerProbe+j])+"\n");
			file.write(repr(B[i*NumberOfCoefficientsPerProbe+j])+"\n");

	print("Saved results to "+StructurePath);
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
	A = np.dot(WeightMatrix.getI(), Colors).T;
	# A = np.dot(pinv(WeightMatrix), ReferenceColors).T;
	B = A.tolist();
	saveCoeffsToFile(B[0], B[1], B[2]);
	
	
# Get the lower corner for trilinear interpolation according to the position
def getSmallestCoordForInterpolation(position, step):
	# print(position);
	x = math.floor((position[0] / step)) * step;
	y = math.floor((position[1] / step)) * step;
	z = math.floor((position[2] / step)) * step;
	
	return [x,y,z]

# Get 
def getInterpolationProbesCoords(position, step):
	smallestCoord = getSmallestCoordForInterpolation(position, step);
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
	print "Entered fillweightmatrix";
	global NumberOfReadSamples, NumberOfProbes, NumberOfCoefficientsPerProbe, WeightMatrix, Colors;
	global index_u, index_v;
	global SurfaceNormalsIJ
	global scene
	step = scene.step;
	sampleIndex = 0;
	# for sampleIndex in xrange(NumberOfReadSamples):
	while (len(WeightMatrix) < NumberOfReadSamples):
		# print("\rRead sample " + repr(sampleIndex+1) + " / " + repr(NumberOfReadSamples)),
		line = GlobalList[:];
		
		# Obtain the Sample Position
		position = SamplePositions[sampleIndex];
		
		# The amounts of probes in each direction and the positions of the first probes in the array
		# dimensions = [21,11,21];
		firstProbePosition = scene.firstProbePosition;
		lastProbePosition = scene.lastProbePosition;
		dimensions = scene.dimensions;

		# Get the positions of the probes to use for interpolation
		coords = getInterpolationProbesCoords(position, step);	

		# Get the importances of each maximal probes for interpolation (upper corners)
		ratios = [(position[0] - coords[0][0]) / step, (position[1] - coords[0][1]) / step, (position[2] - coords[0][2]) / step];

		x = firstProbePosition[0];
		
		count = 0;
		weightSum = 0;

		# Iterate over all the interpolating probes
		for coord in coords:
			# print(("Interpolation probe coord: ", coord))
			badCoord = False;
			
			#Check if calculated interpolation position is outside level
			for i in xrange(3):
				if ((coord[i] > lastProbePosition[i]) or (coord[i] < firstProbePosition[i])):
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
			# ZScene
			# index = int((2*coord[0] + 10) * dimensions[1] * dimensions[2] + \
					# 2 * coord[1] * dimensions[2] + \
					# (2*coord[2] + 10)) * NumberOfCoefficientsPerProbe;

			index = int((coord[0] - firstProbePosition[0]) * dimensions[1] * dimensions[2] + \
						(coord[1] - firstProbePosition[1]) * dimensions[2] + \
						(coord[2] - firstProbePosition[2])) * NumberOfCoefficientsPerProbe;
					
			if (count == 0):
				DFrequencies[index / NumberOfCoefficientsPerProbe] += 1;
			
			# if (DFrequencies[index / NumberOfCoefficientsPerProbe] > 9):
				# break;
			
			# Iterate over all the SH bands
			lm = 0;
			for i in xrange(index, index + NumberOfCoefficientsPerProbe):
				(l,m) = helper.kToLM_new(lm);
				
				# Apply all multipliers for dot product
				# factors = weight * helper.NDotOmegaCoeff(l) * helper.SHxyz_yup((l,m), SampleNormals[sampleIndex]);
				factors = weight * helper.phongCoeffs(l, 1.0) * helper.SHxyz_yup((l,m), SampleNormals[sampleIndex]);
				# print(count, lm, i, factors, weight, helper.NDotOmegaCoeff(l), helper.SHxyz_yup((l,m), SampleNormals[sampleIndex]));

				# print factors
				try:
					line[i] = factors;
				except (IndexError):
					print "Index error!"
					print "at index="+repr(i);
					print "Array length is "+repr(len(line));
					exit();
				lm += 1;
							
			count += 1;

		WeightMatrix.append(line);
		Colors.append(ReferenceColors[sampleIndex]);
		sampleIndex += 1;
	# Get teh matrix from the list of lists
	# WeightMatrix = lil_matrix(WeightMatrix);
	WeightMatrix = np.matrix(WeightMatrix);
	# print(Colors);
	Colors = np.matrix(Colors);

	
	
	
	
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

###########
# Options #
###########
loadRandomProbes = False;
loadFittedProbes = False;
fitCoeffsByInversion = False;
fitCoeffsByMinimizer = False;
showCostFunctions = False;
useRandomWeights = False;
useTriLerpWeights = False;
useUnitWeights = False;
logProbeCoefficients = False;
useExperimentalWeights = False;
useBoxScene = False;

if not (len(sys.argv) > 1):
	print("Possible arguments: ");
	print("fitI: Fit coefficients by matrix pseudo-inversion");
	print("fitM: Fit coefficients by minimizer");
	print("fittedprobes: Use fitted coefficients for the selected parameters");
	print("showcost: Compute and show cost functions");
	print("randomprobes: Load probes with random coefficients");
	print("randomweights: Use random interpolation weights");
	print("trilerp: Use trilinear interpolation weights");
	print("unitweights: Use unitary weights (1 everywhere)");
	print("experimentalweights: Use weights of 1 when i ==j, 0 else");
	print("logcoeffs: Print loaded (not fitted) probe coeffs to the log");
	print("cbox: use cornell box scene");
	print("sponza: use sponza scene");
	sys.exit();

scene = 0;

# reqStep = 1;
reqFirstProbePosition = [0,0,0];
reqDimensions = [7,7,7];
reqStep = 1;
# reqFirstProbePosition = [-10,-1,-7];
# reqDimensions = [21,10,14];

# reqStep = 1;
# reqFirstProbePosition = [-5,0,-5];
# reqDimensions = [11,6,11];

reqFinalProbePos = [sum(x)-1 for x in zip(reqFirstProbePosition, reqDimensions)]
print(reqFinalProbePos)

for arg in sys.argv[1:]:
	if (arg == "fitI"):
		fitCoeffsByInversion = True;
	elif (arg == "fitM"):
		fitCoeffsByMinimizer = True;
	elif (arg == "fittedprobes"):
		loadFittedProbes = True;
	elif (arg == "showcost"):
		showCostFunctions = True;
	elif (arg == "randomprobes"):
		loadRandomProbes = True;
	elif (arg == "randomweights"):
		useRandomWeights = True;
	elif (arg == "trilerp"):
		useTriLerpWeights = True;
	elif (arg == "unitweights"):
		useUnitWeights = True;
	elif (arg == "logcoeffs"):
		logProbeCoefficients = True;
	elif (arg == "experimentalweights"):
		useExperimentalWeights = True;
	elif (arg == "cbox"):
		scene = Scene("zcbox", reqFirstProbePosition, reqFinalProbePos, reqDimensions, 1, 0.01);
	elif (arg == "sponza"):
		scene = Scene("crytek-sponza", reqFirstProbePosition, [10,8,6], reqDimensions, reqStep, 1);
	elif (arg == "zscene"):
		scene = Scene("zscene", reqFirstProbePosition, [5,5,5], reqDimensions, reqStep, 1)
	else:
		print("Unknown argument: "+arg);
		exit();

# Clear log
clearLog();
if (platform.system() == 'Windows'):
	DataFilesPath = 'C:/libraries/g3d/samples/aTest/data-files/';
else:
	DataFilesPath = '/Users/joelpp/Documents/Maitrise/Code/G3D/G3D10/samples/aTest/data-files/';
# print(scene.name);

#########################
# Initialize Parameters #
#########################
NumberOfProbes = scene.dimensions[0] * scene.dimensions[1] * scene.dimensions[2];
NumberOfCoefficientsPerProbe = 9;
NumberOfElements = NumberOfProbes * NumberOfCoefficientsPerProbe;
print(NumberOfElements);

NumberOfReadSamples = 100000;

NumberOfLinesPerSample = 8;

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

######################################
# Load the surface samples from file #
######################################

sampleSetName = "100000_Samples_2";
# sampleSetName = "rgb";

print("Loading samples from " + sampleSetName);

loadSamples(sampleSetName);

###############################
# Load the probe Coefficients #
###############################
probeSetName = "";
if (loadFittedProbes):
	probeSetName = helper.buildStructurePath("Fitted_" + repr(NumberOfReadSamples) + "_Samples", scene.step, scene.dimensions, scene.firstProbePosition);
	loadProbeCoefficients(probeSetName);

elif (loadRandomProbes):
	probeSetName = helper.buildStructurePath("RandomCoefficients", scene.step, scene.dimensions, scene.firstProbePosition);
	loadProbeCoefficients(helper.buildStructurePath("RandomCoefficients", scene.step, scene.dimensions, scene.firstProbePosition));

elif (useBoxScene):
	loadProbeCoefficients("BoxScene_Regular_step1");
else:
	loadProbeCoefficients("RegularGrid_step1");
	
print("Loaded probes from " + probeSetName);

############################################################################################	
# Fill the weight matrix with appropriate interpolation weights / SH functions evaluations #
############################################################################################
if (useRandomWeights):
	WeightMatrix = np.matrix(randomWeights());
	
elif (useUnitWeights):
	for j in xrange(NumberOfReadSamples):
		WeightMatrix.append([]);
		for i in xrange(NumberOfElements):
			weight = 1
			WeightMatrix[j].append(weight);
	WeightMatrix = np.matrix(WeightMatrix);
	
elif (useExperimentalWeights):
	for j in xrange(NumberOfReadSamples):
		WeightMatrix.append([]);
		for i in xrange(NumberOfElements):
			weight = 0;
			if (i == j):
				weight = 1;
			WeightMatrix[j].append(weight);
	WeightMatrix = np.matrix(WeightMatrix);
	
elif (useTriLerpWeights):
	fillWeightMatrix();

print("Filled weight matrix");
# LogPrint("WeightMatrix", WeightMatrix);
# LogPrint("ProbeCoefficientsVector", ProbeCoefficientsVector);

# largest= 0;
# LogPrintList("DFrequencies",DFrequencies);
# for i in xrange(len(DFrequencies)):
	# if (DFrequencies[i] > largest):
		# largest = DFrequencies[i];
		
# print("Largest: " + repr(largest));
# print(WeightMatrix.shape);
# print(ReferenceColors.shape);

# ReferenceColors = np.matrix(ReferenceColors);

# print(ReferenceColors[:,0]);
# result = spsolve(WeightMatrix, ReferenceColors.transpose()[0,:]);

# print(result);




# Fit coefficients to reconstruct reference colors via matrix inversion
if (fitCoeffsByInversion):
	getFittedCoeffs();


if (fitCoeffsByMinimizer):

	COLOR = 0;
	resR = fmin_bfgs(computeRGBs, ProbeCoefficientsVector.T[:,0], maxiter = 1);
	LogPrint("resR",resR);

	COLOR = 1;
	# resG = fmin(computeRGBs, ProbeCoefficientsVector.T[:,1], xtol=1e-10, ftol=1e-10, maxiter = 1e8, maxfun = 1e8);
	resG = fmin_bfgs(computeRGBs, ProbeCoefficientsVector.T[:,1], maxiter = 1);
	LogPrint("resG",resG);

	COLOR = 2;
	# resB = fmin(computeRGBs, ProbeCoefficientsVector.T[:,2], xtol=1e-10, ftol=1e-10, maxiter = 1e8, maxfun = 1e8);
	resB = fmin_bfgs(computeRGBs, ProbeCoefficientsVector.T[:,2], maxiter = 1);
	LogPrint("resB",resB);
	
	saveCoeffsToFile(resR, resG, resB);

if (logProbeCoefficients):
	LogPrint("Probe coeffs", ProbeCoefficientsVector.T);
	
# Compute cost functions with fitted probe coefficients
if (showCostFunctions):
	print("Computing cost functions");
	
	# Compute the cost function in each channel
	COLOR = 0;
	computeRGBs(ProbeCoefficientsVector.T[:,0])
	COLOR = 1
	computeRGBs(ProbeCoefficientsVector.T[:,1])
	COLOR = 2
	computeRGBs(ProbeCoefficientsVector.T[:,2])
	LogPrint("Reconstruction", np.dot(WeightMatrix, ProbeCoefficientsVector.T))
