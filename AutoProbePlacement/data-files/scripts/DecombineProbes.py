import random;
import os;
import helper;

NumberOfReadSamples = 100;
NumberOfCoefficientsPerProbe = 9;

# sceneName = "ZScene"
# firstProbePosition = [-5,0,-5];
# dimensions = [11,6,11];
# structureName = "0"

# sceneName = "crytek-sponza"
# firstProbePosition = [-10,-1,-7];
# dimensions = [21,10,14];
# structureName = "Small_exp";
# firstProbePosition = [1, 4, 0];
# dimensions = [1,1,1];

sceneName = "zcbox"
firstProbePosition = [0,0,0];
dimensions = [7,7,7];
structureName = "Small_exp";

NumberOfProbes = dimensions[0] * dimensions[1] * dimensions[2];
NumberOfElements = NumberOfProbes * NumberOfCoefficientsPerProbe;
print((NumberOfProbes, NumberOfCoefficientsPerProbe))


# StructurePath = "../Scenes/" + sceneName + "/ProbeStructures/" + "Fitted_25000_Samples2";
StructurePath = "../Scenes/" + sceneName + "/ProbeStructures/" + structureName;
# StructurePath = "../Scenes/" + sceneName + "/ProbeStructures/" + "0";
print(StructurePath);


BaseCoeffsFilesPath = "C:/Temp/";
FittedCoeffsFilePath = "";
coeffFiles0 = open(BaseCoeffsFilesPath + "guess0.txt", "w");
coeffFiles1 = open(BaseCoeffsFilesPath + "guess1.txt", "w");
coeffFiles2 = open(BaseCoeffsFilesPath + "guess2.txt", "w");

for probeID in xrange(NumberOfProbes):

	coeffFile = open(StructurePath + "/coefficients/Coefficients_" + repr(probeID)+ ".txt", 'r');
	
	for a in xrange(NumberOfCoefficientsPerProbe):
		a = coeffFile.readline();
		coeffFiles0.write(repr( float(a) ) + "\n")
		coeffFiles1.write(repr( float(coeffFile.readline() ) ) + "\n")
		coeffFiles2.write(repr( float(coeffFile.readline() ) ) + "\n")