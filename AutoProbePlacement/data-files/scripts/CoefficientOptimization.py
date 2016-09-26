# -*- coding: utf-8 -*-
# 
import numpy as np;
from scipy.optimize import minimize;
import helper
<<<<<<< Updated upstream
import platform;
=======

>>>>>>> Stashed changes


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
	samplesFile = open('../temp/samplePositions.txt','r');
	for i in xrange(NumberOfReadSamples):
		ProbeWeights.append([]);
		for j in xrange(NumberOfProbes):
			ProbeWeights[i].append(0);

		for j in xrange(NumberOfLinesPerSample):
			line = samplesFile.readline();
<<<<<<< Updated upstream
			# print("j = "+repr(j)+" , line = "+repr(line));
=======
			print("j = "+repr(j)+" , line = "+repr(line));
>>>>>>> Stashed changes
			if (j == 3): # Normals
				SampleNormals.append(LineToVectorF(line));
			if (j == 4): # Probe Associations
				ProbeAssociations.append(LineToVectorI(line));
			elif (j == 5): # Probe Weights
				Weights = LineToVectorF(line);
				for k in xrange(4):
					ProbeWeights[i][ProbeAssociations[i][k]] = Weights[k];
			elif (j == 6): # Reference RGB
				ReferenceColors.append(LineToVectorF(line));

		samplesFile.readline(); #skip the empty line

def loadProbeCoefficients():
<<<<<<< Updated upstream

	for i in xrange(NumberOfCoefficientsPerProbe):
		ProbeCoefficientsBlock[0].append([]);
		ProbeCoefficientsBlock[1].append([]);
		ProbeCoefficientsBlock[2].append([]);

	for i in xrange(NumberOfProbes):
		# CoefficientFile = open(DataFilesPath+'ProbeStructures/test/coefficients/Coefficients_'+repr(i)+'.txt');
=======
	for i in xrange(NumberOfProbes):
		# CoefficientFile = open(DataFilesPath+'ProbeStructures/test/coefficients/Coefficients_'+repr(i)+'.txt');
		# ProbeCoefficients.append([]);
>>>>>>> Stashed changes
		# for j in xrange(NumberOfCoefficientsPerProbe):
		# 	ProbeCoefficients[i].append([]);
		# 	for k in xrange(3):
		# 		ProbeCoefficients[i][j].append(float(CoefficientFile.readline()));
		CoefficientFile = open(DataFilesPath+'ProbeStructures/test/coefficients/Coefficients_'+repr(i)+'.txt');
<<<<<<< Updated upstream
		ProbeCoefficients[0].append([]);
		ProbeCoefficients[1].append([]);
		ProbeCoefficients[2].append([]);

		for j in xrange(NumberOfCoefficientsPerProbe):
			# ProbeCoefficients[0][i].append(float(CoefficientFile.readline()));
			# ProbeCoefficients[1][i].append(float(CoefficientFile.readline()));
			# ProbeCoefficients[2][i].append(float(CoefficientFile.readline()));
			ProbeCoefficientsBlock[0][j].append(float(CoefficientFile.readline()));
			ProbeCoefficientsBlock[1][j].append(float(CoefficientFile.readline()));
			ProbeCoefficientsBlock[2][j].append(float(CoefficientFile.readline()));

	helper.show("ProbeCoefficientsBlock", np.bmat(ProbeCoefficientsBlock));
	# for i in xrange(3):
		# ProbeCoefficients[i] = np.matrix(ProbeCoefficients[i]);
		
def init():
	loadSamples();
	loadProbeCoefficients();
	BuildSHMatrix();


def getAlYlmVector(normal):
	toReturn = [];
	for i in xrange(NumberOfCoefficientsPerProbe):
		(l,m) = helper.kToLM_new(i);
		toReturn.append(helper.NDotOmegaCoeff(l) * helper.SHxyz_yup((l,m), normal));

	return np.array(toReturn);
	
def BuildSHMatrix():
	for i in xrange(NumberOfReadSamples):
		SHMatrix.append(getAlYlmVector(SampleNormals[i]));
	# print(SHMatrix);
	global SHMatrix;
	SHMatrix = np.matrix(SHMatrix);
	
def debugPrinting():
	# print(np.array(ReferenceColors));
	printSeparator();
	A = ProbeCoefficients[0].T;
=======
		ProbeCoefficientsR.append([]);
		ProbeCoefficientsG.append([]);
		ProbeCoefficientsB.append([]);

		for j in xrange(NumberOfCoefficientsPerProbe):
			ProbeCoefficientsR[i].append(float(CoefficientFile.readline()));
			ProbeCoefficientsG[i].append(float(CoefficientFile.readline()));
			ProbeCoefficientsB[i].append(float(CoefficientFile.readline()));


def init():
	loadSamples();
	loadProbeCoefficients();


def cost(x):
	# return la somme des diff au carré
	# x vecteur d'entrée : les coefficients des probes
	# 
	calculated = [];
	# print(x);

	for i in xrange(NumberOfReadSamples):      
		calculated.append([0,0,0]);

	for i in xrange(NumberOfReadSamples):
		for j in xrange(4):
			ProbeToUse = ProbeAssociations[i][j];
			weight = ProbeWeights[i][j];
			
			# print("ProbeToUse = "+repr(ProbeToUse));
			# print("weight = "+repr(weight));
			for k in xrange(NumberOfSHBands):
				# print(x[ProbeToUse][k]);
				calculated[i][0] += x[ProbeToUse][k][0] * weight;
				calculated[i][1] += x[ProbeToUse][k][1] * weight;
				calculated[i][2] += x[ProbeToUse][k][2] * weight;

	# print(np.array(calculated));
	# print(np.array(ReferenceColors));
	# print(np.array(calculated) - np.array(ReferenceColors));
	# print( np.sum(np.array(calculated) - np.array(ReferenceColors), axis=1)**2);
	return sum( np.sum(np.array(calculated) - np.array(ReferenceColors), axis=1)**2);

def debugPrinting():
	# print(np.array(ReferenceColors));
	printSeparator();
	A = np.matrix(ProbeCoefficientsR).T;
>>>>>>> Stashed changes
	B = np.array(ProbeWeights);
	print(A);
	printSeparator();
	# print(ProbeAssociations);
	print(B);
	printSeparator();
	print(np.dot(A,B.T));
	printSeparator();
	print(SampleNormals);
	printSeparator();
	# print(getAlYlmVector(SampleNormals[0]));
	printSeparator();
	print(np.dot(A,B.T).T.dot(getAlYlmVector(SampleNormals[0])));
<<<<<<< Updated upstream
	print(np.dot(A,B.T).T.dot(getAlYlmVector(SampleNormals[1])));

	print(np.dot(A,B.T).T.dot(SHMatrix.T));
	
def cost(x):
	# return la somme des diff au carré
	# x vecteur d'entrée : les coefficients des probes
	global PassIndex;
	
	A = x.T;
	B = np.array(ProbeWeights);
	PassIndex = PassIndex + 1;
	InterpolatedCoefficients =  np.dot(A,B.T);
	result = np.diagonal(InterpolatedCoefficients.T.dot(SHMatrix.T));
	# helper.show("PassIndex", PassIndex);
	helper.show("ProbeCoefficients (x.T)", x.T);
	helper.show("ProbeWeights", B);
	helper.show("InterpolatedCoefficients", InterpolatedCoefficients);

	helper.show("SHMatrix.T", SHMatrix.T);

	helper.show("ReferenceColors", ReferenceColors[:,0]);
	helper.show("final:", result);
	return ;
	

	

# Program Start	
NumberOfProbes = 9;
NumberOfCoefficientsPerProbe = 9;
NumberOfReadSamples = 2;
=======
def getAlYlmVector(normal):
	toReturn = [];
	for i in xrange(NumberOfCoefficientsPerProbe):
		(l,m) = helper.kToLM_new(i);
		print (l, m);
		print helper.NDotOmegaCoeff(l);
		toReturn.append(helper.NDotOmegaCoeff(l) * helper.SHxyz_yup((l,m), normal));

	return np.array(toReturn);
# Program Start	
NumberOfProbes = 9;
NumberOfCoefficientsPerProbe = 9;
NumberOfReadSamples = 5;
>>>>>>> Stashed changes
NumberOfLinesPerSample = 7;
NumberOfSHBands = 9;
SHCoeffs = [];

<<<<<<< Updated upstream
if (platform.system() == 'Windows'):
	DataFilesPath = 'C:/libraries/g3d/samples/aTest/data-files/';
else:
	DataFilesPath = '/Users/joelpp/Documents/Maitrise/Code/G3D/G3D10/samples/aTest/data-files/';
ReferenceColors = [];
ProbeCoefficients = [];
ProbeCoefficientsR = [];
ProbeCoefficientsG = [];
ProbeCoefficientsB = [];
ProbeCoefficients = [[],[],[]];
ProbeCoefficientsVector = [];
ProbeCoefficientsBlock = [[],[],[]];
=======
DataFilesPath = '/Users/joelpp/Documents/Maitrise/Code/G3D/G3D10/samples/aTest/data-files/';
ReferenceColors = [];
ProbeCoefficients = [];

ProbeCoefficientsR = [];
ProbeCoefficientsG = [];
ProbeCoefficientsB = [];
>>>>>>> Stashed changes
ProbeAssociations = [];
ProbeWeights = [];
SampleNormals = [];
DummyRGBColors = [];
<<<<<<< Updated upstream
SHMatrix = [];
PassIndex = 0;

init();
# debugPrinting();
printSeparator();
printSeparator();

R = cost(ProbeCoefficients[0]);
# G = cost(ProbeCoefficients[1]);
# B = cost(ProbeCoefficients[2]);
# print R,G,B
# result = np.matrix([R,G,B]);

print R;
=======

init();
debugPrinting();
# print(cost(ProbeCoefficients));

>>>>>>> Stashed changes
# TestArray = np.array([[1,1,1], [0.2,0.2,0.2]]);
# printSeparator();

# print(TestArray);
# printSeparator();
# print(TestArray - ReferenceColors);
# printSeparator();
# print((TestArray - ReferenceColors)**2);

# x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2]);
<<<<<<< Updated upstream
res = minimize(cost, ProbeCoefficients[0], method='nelder-mead',options={'xtol': 1e-8, 'disp': True});
=======
# res = minimize(cost, ProbeCoefficients, method='nelder-mead',options={'xtol': 1e-8, 'disp': True});
>>>>>>> Stashed changes
# print(res);