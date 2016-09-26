# -*- coding: utf-8 -*-
# 
from __future__ import division
import winsound
import numpy as np;
from scipy.optimize import minimize;
import helper
import platform;
np.set_printoptions(precision=5)
np.set_printoptions(suppress=True)
from scipy import misc
import sys
import math
import multiprocessing
LogPath = "../temp/SurfaceTextureGenerationLog.txt";
import os


def clearLog():
	writeStringToFile(LogPath, "", 'w');

def Lolprint(MatrixName, Matrix):
	try:
		writeStringToFile(LogPath, MatrixName+" \n" + assembleMatrixString(Matrix.tolist()));
	except:
		try:
			writeStringToFile(LogPath, MatrixName+" \n" + assembleMatrixString(Matrix));
		except:
			writeStringToFile(LogPath, MatrixName+" \n" + repr(Matrix) + '\n');

def writeStringToFile(Path, String, option='a'):
	file = open(Path, option);
	file.write(String+"\n");
	file.close();
	
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

def printSeparator():
	print("---------------------------");


def rosen(x):
	# """The Rosenbrock function""";
	return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0);

def LineToVectorF(line):
	vector = [float(x) for x in line.split()]
	# print(vector);
	return (vector);

def LineToVectorI(line):
	vector = [int(x) for x in line.split()]
	# print(vector);
	return (vector);

def AppendToArray(array, vector):
	array.append(vector);

def loadSamples():
	global ReferenceColors;
	global Values;
	if (useBoxScene):
		samplesFile = open('../temp/surfaceSamples/128x128/surfaceSamples_'+repr(sampleID)+'.txt','r');
	else:
		samplesFile = open('../temp/zScene/surfaceSamples/surfaceSamples_'+repr(sampleID)+'.txt','r');

	for i in xrange(NumberOfReadSamples):
		print "\rLoading Sample " + repr(i) + " / " + repr(NumberOfReadSamples),;

		for j in xrange(NumberOfLinesPerSample):
			line = samplesFile.readline();
			# print("j = "+repr(j)+" , line = "+repr(line));
			if (j == 1): # Normals
				a = LineToVectorF(line);
				SamplePositions.append(a);
				
				SampleNormals.append(SurfaceNormals[sampleID]);
			elif (j == 2): # Probe Associations
				# ProbeAssociations.append(LineToVectorI(line));
				a = LineToVectorF(line);
				ReferenceColors.append(a);
			# elif (j == 3): # Probe Weights
				# Weights = LineToVectorF(line);
				# for k in xrange(4):
					# ProbeWeights[i][ProbeAssociations[i][k]] = Weights[k];
			# elif j == 4:
				# Values.append(float(line));

		samplesFile.readline(); #skip the empty line

	print "\n";
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

	
def getSmallestCoordForInterpolation(position):
	x = math.floor((position[0] / step)) * step;
	y = math.floor((position[1] / step)) * step;
	z = math.floor((position[2] / step)) * step;
	return [x,y,z]
	
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

def fillWAYMatrix3(sampleIndex):
	global NumberOfReadSamples, NumberOfProbes, NumberOfCoefficientsPerProbe, WAYMatrix;
	line = GlobalList[:];
	global index_u, index_v;
	global SurfaceNormalsIJ
	# Start out by filling the line with 0s

	
	# if (ReferenceColors[sampleIndex] == [0,0,0]):
		# return line;
	# print(("total size of line", len(line)));
	
	# Obtain the Sample Position
	position = SamplePositions[sampleIndex];
	# position = [-4.09411, 2.30075, -3.9282]
	# position = [-3.78114, 3.46798, -3.6747]
	
	# The amounts of probes in each direction and the positions of the first probes in the array
	# repetitions = [21,11,21];
	
	if (useBoxScene):
		firstProbePosition = [-2,0,-2];
		lastProbePosition = [2,2,2];
		repetitions = [5,3,5];

	else:
		firstProbePosition = [-5,0,-5];
		lastProbePosition = [5,5,5];
		repetitions = [11,6,11];
	

	# Get the positions of the probes to use for interpolation
	coords = getInterpolationProbesCoords(position);	
	if (coords[0] not in SeenMinCoords):
		SeenMinCoords.append(coords[0]);
	# return SeenMinCoords.index(coords[0]);
	# Get the importances of each maximal probes for interpolation (upper corners)
	ratios = [(position[0] - coords[0][0]) / step, (position[1] - coords[0][1]) / step, (position[2] - coords[0][2]) / step];
	# return ratios
	

	x = firstProbePosition[0];
	count = 0;
	weightSum = 0;
	if (sampleIndex == 883):
		Lolprint("position",position)
		Lolprint("ratios",ratios)
		Lolprint("coords",coords)

	# Iterate over all the interpolating probes
	for coord in coords:
		# print(("Interpolation probe coord: ", coord))
		badCoord = False;
		for i in xrange(3):
			if (coord[i] > lastProbePosition[i]):
				badCoord = True;
				break;
		if (badCoord == True):
			count += 1;
			continue;
			
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
		if (useBoxScene):
			index = int((coord[0] + 2) * repetitions[1] * repetitions[2] + \
						coord[1] * repetitions[2] + \
						(coord[2] + 2)) * NumberOfCoefficientsPerProbe;
		else:
			index = int((coord[0] + 5) * repetitions[1] * repetitions[2] + \
						coord[1] * repetitions[2] + \
						(coord[2] + 5)) * NumberOfCoefficientsPerProbe;
						
		# ZScene
		# index = int((2*coord[0] + 10) * repetitions[1] * repetitions[2] + \
				# 2 * coord[1] * repetitions[2] + \
				# (2*coord[2] + 10)) * NumberOfCoefficientsPerProbe;

		

		lm = 0;
		
		# Iterate over all the SH bands
		for i in xrange(index, index + NumberOfCoefficientsPerProbe):
			(l,m) = helper.kToLM_new(lm);
			
			# Apply all multipliers for dot product

			try:
				line[i] = weight * helper.NDotOmegaCoeff(l) * helper.SHxyz_yup((l,m), SurfaceNormals[sampleID]);
			except (IndexError):
				print "Something wrong with your indexing!"
				print "SampleID = "+repr(sampleID);
				print "i = "+repr(i);
				print "coord = "+repr(coord);
				exit();

			lm += 1;
			
		# print "\n";
		
		count += 1;

	# Return the line with all proper multiplicators inserted
	return line;


def saveCoeffsToFile():
	file = open(DataFilesPath+"ProbeStructures/computed/GeneratedCoeffFromSamples.txt", 'w');

	for i in xrange(NumberOfCoefficientsPerProbe):
		# file = open(DataFilesPath+"ProbeStructures/computed/coefficients/Coefficients_"+repr(i)+".txt", 'w');
		print(ComputedCoeffs[0][i]);
		file.write(repr(ComputedCoeffs[0][i])+"\n");
		file.write(repr(ComputedCoeffs[1][i])+"\n");
		file.write(repr(ComputedCoeffs[2][i])+"\n");


		
def computeColors():
	# Initialize samples matrix to 0
	lines = [];
	
	# For each sample (i.e. in this case, texel on the texture being generated
	for i in xrange(NumberOfReadSamples):
	
		print "\rSurface "+repr(sampleID)+" : Computing sample " + repr(i) + " / " + repr(NumberOfReadSamples),;
	
		# Fill a line corresponding to a texel with interpolation weight * geometric term * SH basis function evaluation in direction of normal
		
		
		line = fillWAYMatrix3(i);
		
		# Convert the line to its own matrix so we can do dot product
		line = np.matrix(line);
		# Lolprint(repr(i), line.T);
		
		# dot it with the ProbeCoefficients, convert to a list and add the resulting color
		# to texel colors matrix
		result = np.dot(line, ProbeCoefficientsVector.T).tolist()[0];
		# result = [line/50.0,line/50.0,line/50.0];
		# result = line;
		# Lolprint("result"+repr(i), result);
		lines.append(result);
	print "\n"

	# For testing - directly compute the interpolated coefficients and save them to a file
	# for i in xrange(NumberOfCoefficientsPerProbe):
		# ComputedCoeffs[0].append(0);
		# ComputedCoeffs[1].append(0);
		# ComputedCoeffs[2].append(0);

		# for j in xrange(NumberOfProbes):
			# index = i*NumberOfProbes+j;
			# ComputedCoeffs[0][i] += lines[0][i+NumberOfCoefficientsPerProbe*j] * ProbeCoefficientsVector[0][i+NumberOfCoefficientsPerProbe*j];
			# ComputedCoeffs[1][i] += lines[0][i+NumberOfCoefficientsPerProbe*j] * ProbeCoefficientsVector[1][i+NumberOfCoefficientsPerProbe*j];
			# ComputedCoeffs[2][i] += lines[0][i+NumberOfCoefficientsPerProbe*j] * ProbeCoefficientsVector[2][i+NumberOfCoefficientsPerProbe*j];
			
	return lines;
	
# General purpose debug printing
def debugPrinting():
	Lolprint("ProbeCoefficientsVector", ProbeCoefficientsVector);
	Lolprint("ProbeCoefficientsVector[0]", ProbeCoefficientsVector[0]);
	Lolprint("ProbeCoefficientsVector[0][0]", ProbeCoefficientsVector[0][0]);
	Lolprint("ProbeCoefficientsVector[0][1]", ProbeCoefficientsVector[0][1]);
	Lolprint("ProbeCoefficientsVector[1][0]", ProbeCoefficientsVector[1][0]);
	Lolprint("ProbeCoefficientsVector.T", ProbeCoefficientsVector.T);

	pass
	
#################
# Program Start	#
#################

###########
# Options #
###########

useBoxScene = False;

if not (len(sys.argv) > 1):
	print("useboxscene: use small box scene");
	sys.exit();

for arg in sys.argv[1:]:
	if (arg == "useboxscene"):
		useBoxScene = True;
	elif (arg == "usezscene"):
		useBoxScene = False;
	else:
		print("Unknown argument: "+arg);
		exit();

# Sample set data
if (useBoxScene):
	SurfaceNormals = [[1,0,0],[0,0,1],[-1,0,0],[0,0,-1],[0,-1,0],[0,1,0]]; #boxscene
	NumberOfProbes = 74;
	textureWidth = 128;
	textureHeight = 128;
else:
	SurfaceNormals = [[-1,0,0],[0,0,1],[1,0,0],[0,0,-1],[0,1,0],[0,-1,0],[1,0,0],[0,0,1],[-1,0,0],[0,0,-1]]; #zscene
	NumberOfProbes = 726;
	textureWidth = 128;
	textureHeight = 64;


# SurfaceNormals = [[],[],[]]
SurfaceNormalsIJ = [[128,0],[128, 384]]; #unflipped
albedos = [[1,1,1],[0,0,1],[1,1,1],[0,0,1],[1,1,1],[1,1,1],[1,0,0],[1,0,0],[0,1,0],[0,1,0]];
flippedalbedos = [[1,1,1],[1,0,0],[1,1,1],[1,0,0],[1,1,1],[1,1,1],[0,0,1],[0,0,1],[0,1,0],[0,1,0]];

Freq = 2500 # Set Frequency To 2500 Hertz
Dur = 200 # Set Duration To 1000 ms == 1 second

list = [];

ProbeCoefficientsVector = [[],[],[]];
if (platform.system() == 'Windows'):
	DataFilesPath = 'C:/libraries/g3d/samples/aTest/data-files/';
else:
	DataFilesPath = '/Users/joelpp/Documents/Maitrise/Code/G3D/G3D10/samples/aTest/data-files/';

	clearLog();

# Probe data
NumberOfCoefficientsPerProbe = 9;	
step = 1;


NumberOfReadSamples = textureWidth * textureHeight;
NumberOfLinesPerSample = 4;
NumberOfSurfaceSamples = 90000;

# Chck if these textures have been calculated
if (useBoxScene):	
	structureID = "BoxScene_computed_"+repr(NumberOfSurfaceSamples)+"samples_"+repr(NumberOfCoefficientsPerProbe)+"coeffs/";
else:
	structureID = "computed_"+repr(NumberOfSurfaceSamples)+"samples_"+repr(NumberOfCoefficientsPerProbe)+"coeffs/";
savePath = "../temp/surfaceRenders/"+structureID;
if (os.path.isdir(savePath)):
	choice = raw_input("This configuration exists - are you sure you want to overwrite? (y/n)");
	if not (choice == 'y'):
		exit();
else:
	os.makedirs(savePath);
	os.makedirs(savePath+"Reconstruction_optimized/");
	os.makedirs(savePath+"Reference/");
	os.makedirs(savePath+"X/");
	os.makedirs(savePath+"Y/");
	os.makedirs(savePath+"Z/");

# Debugging arrays
SeenMinCoords = [];
ComputedCoeffs = [[],[],[]];

# Load and store the probe set
loadProbeCoefficients(structureID);
GlobalList = [];
for i in xrange(NumberOfProbes * NumberOfCoefficientsPerProbe):
	GlobalList.append(0);
	
# Iterate on each surface in the scene
if (useBoxScene):
	numSurfaces = 6;
else:
	numSurfaces = 10;
for sampleID in xrange(numSurfaces):
	SeenMinCoords = [];

	# Surface specific data
	SamplePositions = [];
	SampleNormals = [];
	ReferenceColors = [];

	# Load and store the samples
	loadSamples();
	printSeparator();
	
	# Print debug stuff if you feel like it
	# debugPrinting();

	# Compute the colors array for the surface
	c = computeColors();
	Lolprint("c", c);
	
	# Save two versions to file for comparison, one multiplied by albedo color and one not
	RGB = [];
	RGB2 = [];
	
	# Also save a bunch of debugging images
	RGBPositions = [];
	RGBX = [];
	RGBY = [];
	RGBZ = [];
	RGBReference = [];
	print "finished computing colors";
	
	for heightIndex in xrange(textureHeight):
		RGB.append([]);
		RGB2.append([]);
		RGBPositions.append([]);
		RGBReference.append([]);
		RGBX.append([]);
		RGBY.append([]);
		RGBZ.append([]);
	for heightIndex in xrange(textureHeight):

		for widthIndex in xrange(textureWidth):
			# print(c[heightIndex*textureWidth+widthIndex]);
			index = heightIndex*textureWidth+widthIndex;
			height = textureHeight - heightIndex - 1;

			RGB[height].append([]);
			RGB2[height].append([]);
			

			#Debug stuff
			# RGBPositions[height].append([(SamplePositions[index][0]+5.0)/10.0, \
								  # SamplePositions[index][1]/5.0, \
								  # (SamplePositions[index][2]+5.0)/10.0]);
			# RGBX[height].append((SamplePositions[index][0] + 5) / 5.0);
			# RGBY[height].append((SamplePositions[index][1]) / 5.0);
			# RGBZ[height].append((SamplePositions[index][2] + 5) / 5.0);
			RGBReference[height].append(ReferenceColors[index]);
			
			#Good stuff
			RGB[height][widthIndex] = c[index];
			# RGB2[height][widthIndex] = [a*b/3.141592654 for a,b in zip(c[index],albedos[sampleID])];
			# RGB[heightIndex][widthIndex] = ReferenceColors[heightIndex*textureWidth+widthIndex];
			# RGB2[heightIndex][widthIndex] = [a*b/3.141592654 for a,b in zip(ReferenceColors[heightIndex*textureWidth+widthIndex],albedos[sampleID])];

	# misc.imsave("../temp/surfaceRenders/c"+repr(sampleID)+".png", c);


	im = misc.toimage(RGB, cmin=0, cmax = 1);
	im.save(savePath+"surfaceRender_"+repr(sampleID)+".png");
	im.save(savePath+"Reconstruction_optimized/surfaceRender_"+repr(sampleID)+".png");
	print "Saved reconstruction to "+savePath+"surfaceRender_"+repr(sampleID)+".png";
	
	# misc.imsave("../temp/surfaceRenders/surfaceRender_"+repr(sampleID)+".png", RGB);
	# misc.imsave(savePath+"surfaceRenderAlbedo_"+repr(sampleID)+".png", RGB2);
	# misc.imsave(savePath+"surfaceRenderPosition_"+repr(sampleID)+".png", RGBPositions);
	
	misc.imsave(savePath+"surfaceRenderReference_"+repr(sampleID)+".png", RGBReference);
	misc.imsave(savePath+"Reference/surfaceRenderReference_"+repr(sampleID)+".png", RGBReference);

	# misc.imsave(savePath+"X/surfaceRenderX_"+repr(sampleID)+".png", RGBX);
	# misc.imsave(savePath+"surfaceRenderX_"+repr(sampleID)+".png", RGBX);
	
	# misc.imsave(savePath+"Y/surfaceRenderY_"+repr(sampleID)+".png", RGBY);
	# misc.imsave(savePath+"surfaceRenderY_"+repr(sampleID)+".png", RGBY);

	# misc.imsave(savePath+"Z/surfaceRenderZ_"+repr(sampleID)+".png", RGBZ);
	# misc.imsave(savePath+"surfaceRenderZ_"+repr(sampleID)+".png", RGBZ);

# winsound.Beep(Freq,Dur)