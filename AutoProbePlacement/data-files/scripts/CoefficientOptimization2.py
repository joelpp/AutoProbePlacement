# -*- coding: utf-8 -*-
# 
import numpy as np;
from scipy.optimize import minimize;
import helper
import platform;
np.set_printoptions(precision=5)
np.set_printoptions(suppress=True)

LogPath = "../temp/PythonLog.txt";

global COLOR;
COLOR = -1;

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
	samplesFile = open('../temp/samplePositions.txt','r');
	for i in xrange(NumberOfReadSamples):
		ProbeWeights.append([]);
		for j in xrange(NumberOfProbes):
			ProbeWeights[i].append(0);

		for j in xrange(NumberOfLinesPerSample):
			line = samplesFile.readline();
			# print("j = "+repr(j)+" , line = "+repr(line));
			if (j == 0): # Surface IDs
				SampleAlbedos.append(albedos[LineToVectorI(line)[0]]);
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
	ReferenceColors = np.matrix(ReferenceColors);
	LogPrint("ReferenceColors", ReferenceColors);
	# LogPrint("ProbeWeights", ProbeWeights);
	
def loadProbeCoefficients():
	global ProbeCoefficientsVector;
	for i in xrange(NumberOfCoefficientsPerProbe):
		for j in xrange(NumberOfProbes):
			ProbeCoefficientsVector[0].append(0);
			ProbeCoefficientsVector[1].append(0);
			ProbeCoefficientsVector[2].append(0);

	for i in xrange(NumberOfProbes):
		CoefficientFile = open(DataFilesPath+'ProbeStructures/currenter_morecoeffs/coefficients/Coefficients_'+repr(i)+'.txt');

		for j in xrange(NumberOfCoefficientsPerProbe):

			ProbeCoefficientsVector[0][i*NumberOfProbes+j] = float(CoefficientFile.readline());
			ProbeCoefficientsVector[1][i*NumberOfProbes+j] = float(CoefficientFile.readline());
			ProbeCoefficientsVector[2][i*NumberOfProbes+j] = float(CoefficientFile.readline());

	ProbeCoefficientsVector = np.array(ProbeCoefficientsVector);
	LogPrint("ProbeCoefficientsVector.T", ProbeCoefficientsVector.T);


def init():
	loadSamples();
	loadProbeCoefficients();
	fillWAYMatrix();

def getAlYlmVector(normal):
	toReturn = [];
	for i in xrange(NumberOfCoefficientsPerProbe):
		(l,m) = helper.kToLM_new(i);
		toReturn.append(helper.NDotOmegaCoeff(l) * helper.SHxyz_yup((l,m), normal));

	return np.array(toReturn);
	

	
def clearLog():
	writeStringToFile(LogPath, "", 'w');

def LogPrint(MatrixName, Matrix):
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

	
def fillWAYMatrix():
	print("Entered fillWAYMatrix()");
	global WAYMatrix;
	for i in xrange(NumberOfReadSamples):
		WAYMatrix.append([]);
			
		for j in xrange(NumberOfProbes):
			for k in xrange(NumberOfCoefficientsPerProbe):
				(l,m) = helper.kToLM_new(k);
				# print((i,j,(l,m),ProbeWeights[i][j], helper.NDotOmegaCoeff(l), helper.SHxyz_yup((l,m), SampleNormals[i])))
				WAYMatrix[i].append(ProbeWeights[i][j] * helper.NDotOmegaCoeff(l) * helper.SHxyz_yup((l,m), SampleNormals[i]));
	

	WAYMatrix = np.matrix(WAYMatrix);

	LogPrint("WAYMatrix",WAYMatrix);
	LogPrint("WAYMatrix.getI()", np.linalg.pinv(WAYMatrix));

def saveCoeffsToFile(R, G, B):
	for i in xrange(NumberOfProbes):
		file = open(DataFilesPath+"ProbeStructures/computed/coefficients/Coefficients_"+repr(i)+".txt", 'w');
		for j in xrange(NumberOfCoefficientsPerProbe):
			file.write(repr(R[i*NumberOfProbes+j])+"\n");
			file.write(repr(G[i*NumberOfProbes+j])+"\n");
			file.write(repr(B[i*NumberOfProbes+j])+"\n");

# def computeRGBs():
	# # helper.show("WAYMatrix", WAYMatrix.shape);
	# # helper.show("ReferenceColors[:,0].shape", ReferenceColors[:,0].shape);
	# # helper.show("WAYMatrix", WAYMatrix);
	# helper.show("ReferenceColors", ReferenceColors);
	# A = np.dot(WAYMatrix, ProbeCoefficientsVector.T);
	# print(A);
	# helper.show("diff", A - ReferenceColors);
	# helper.show("diff**2", np.square(A - ReferenceColors));
	# helper.show("sum(diff**2)", np.sum(np.square(A - ReferenceColors)));


	
def computeRGBs(x):
	global PassIndex
	global PRINTLOG
	global COLOR

	A = np.dot(WAYMatrix, x);
	B = A.T*SampleAlbedos[:][COLOR] - ReferenceColors[:,COLOR];

	C = np.square(B);
	D = np.sum(C);
	
	if (PRINTLOG == 1):
		LogPrint("np.dot(WAYMatrix, x).T", A.T);
	# if (PassIndex % 1000 == 0):
		# LogPrint("PassIndex",PassIndex);
		# LogPrint("x",x);
		# LogPrint("A",A);
		# LogPrint("B",B);
		# LogPrint("C",C);	
		# LogPrint("D",D);
	
	PassIndex += 1;
	
	return D;
	

def debugPrinting():
	global A;

	LogPrint("ReferenceColors",ReferenceColors[:,0]);
	LogPrint("np.dot(WAYMatrix, ProbeCoefficientsVector.T[:,0]).T", np.dot(WAYMatrix, ProbeCoefficientsVector.T[:,0]).T);
	LogPrint("np.dot(WAYMatrix, ProbeCoefficientsVector.T).T - ReferenceColors", np.dot(WAYMatrix, ProbeCoefficientsVector.T[:,0]).T - ReferenceColors[:,0]);
	LogPrint("np.sum(np.square(np.dot(WAYMatrix, ProbeCoefficientsVector.T[:,0]).T - ReferenceColors[:,0]))", np.sum(np.square(np.dot(WAYMatrix, ProbeCoefficientsVector.T[:,0]).T - ReferenceColors[:,0])));
			
			
			
def getFittedCoeffs():

	global A;
	A = np.dot(WAYMatrix.getI(), ReferenceColors).T;
	B = A.tolist();
	LogPrint("Computed Coefficients from inversion (transposed)",A.T);
	for i in xrange(NumberOfProbes):
		file = open(DataFilesPath+"ProbeStructures/computed/coefficients/Coefficients_"+repr(i)+".txt", 'w');
		for j in xrange(NumberOfCoefficientsPerProbe):
			file.write(repr(B[0][i*9+j])+"\n");
			file.write(repr(B[1][i*9+j])+"\n");
			file.write(repr(B[2][i*9+j])+"\n");
			
# Program Start	

clearLog();
if (platform.system() == 'Windows'):
	DataFilesPath = 'C:/libraries/g3d/samples/aTest/data-files/';
else:
	DataFilesPath = '/Users/joelpp/Documents/Maitrise/Code/G3D/G3D10/samples/aTest/data-files/';

NumberOfProbes = 9;
NumberOfCoefficientsPerProbe = 9;
NumberOfReadSamples = 81;
NumberOfLinesPerSample = 7;
NumberOfSHBands = 9;
albedos = [[1,1,1],[0,0,1],[1,1,1],[0,0,1],[1,1,1],[1,1,1],[1,0,0],[1,0,0],[0,1,0],[0,1,0]];

WAYMatrix = [];

ReferenceColors = [];

ProbeCoefficientsVector = [[],[],[]];
ProbeAssociations = [];
ProbeWeights = [];
SampleNormals = [];
SampleAlbedos = [];

PassIndex = 0;
PRINTLOG = 0;
init();
# printSeparator();
# printSeparator();
# debugPrinting();
getFittedCoeffs();


# LogPrint("SampleAlbedos",SampleAlbedos);

# LogPrint("computeRGBs(ProbeCoefficientsVector.T[:,0])",computeRGBs(ProbeCoefficientsVector.T[:,0]));
# LogPrint("computeRGBs(ProbeCoefficientsVector.T[:,1])",computeRGBs(ProbeCoefficientsVector.T[:,1]));
# LogPrint("computeRGBs(ProbeCoefficientsVector.T[:,2])",computeRGBs(ProbeCoefficientsVector.T[:,2]));

# resR = minimize(computeRGBs, ProbeCoefficientsVector.T[:,0], method='nelder-mead',options={'xtol': 1e-8, 'disp': True, 'maxfun':92233720368547});
# COLOR = 1;
# resG = minimize(computeRGBs, ProbeCoefficientsVector.T[:,1], method='nelder-mead',options={'xtol': 1e-8, 'disp': True, 'maxfun':92233720368547});
# COLOR = 2;
# resB = minimize(computeRGBs, ProbeCoefficientsVector.T[:,2], method='nelder-mead',options={'xtol': 1e-8, 'disp': True, 'maxfun':92233720368547});

# LogPrint("resR",resR);
# LogPrint("resG",resG);
# LogPrint("resB",resB);

# # LogPrint("WAYMatrix.dot(OptimizedCoeffs.T) - ReferenceColors", np.dot(WAYMatrix, resR.x.T)[:,0] - ReferenceColors[:,0])
# # PassIndex = 0;
# # PassIndex = 0;
# # print(resR.x);
# PRINTLOG = 1;
# LogPrint("computeRGBs(ProbeCoefficientsVector.T[:,0])",computeRGBs(resR.x));
# LogPrint("computeRGBs(ProbeCoefficientsVector.T[:,0])",computeRGBs(resG.x));
# LogPrint("computeRGBs(ProbeCoefficientsVector.T[:,0])",computeRGBs(resB.x));

# saveCoeffsToFile(resR.x, resB.x, resG.x);
