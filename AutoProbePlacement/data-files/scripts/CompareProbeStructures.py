import random;
import os;
import helper;

NumberOfCoefficientsPerProbe = 9;

# sceneName = "ZScene"
# firstProbePosition = [-5,0,-5];
# dimensions = [11,6,11];
# structureName = "0"

# sceneName = "crytek-sponza"
# firstProbePosition = [-10,-1,-7];
# dimensions = [21,10,14];
# structure0Name = "Small_exp";
# structure1Name = "Fitted_5000_Samples0";
# firstProbePosition = [1, 4, 0];
# dimensions = [1,1,1];

sceneName = "zcbox"
firstProbePosition = [0,0,0];
dimensions = [7,7,7];
structure0Name = "Small_exp2";
# structure1Name = "Fitted_100000_Samples4";
structure1Name = "Fitted_150000_Samples_step_1_dimensions_7_7_7_firstProbePosition_0_0_0";

NumberOfProbes = dimensions[0] * dimensions[1] * dimensions[2];
NumberOfElements = NumberOfProbes * NumberOfCoefficientsPerProbe;
print((NumberOfProbes, NumberOfCoefficientsPerProbe))


# StructurePath = "../Scenes/" + sceneName + "/ProbeStructures/" + "Fitted_25000_Samples2";
Structure0Path = "../Scenes/" + sceneName + "/ProbeStructures/" + structure0Name;
Structure1Path = "../Scenes/" + sceneName + "/ProbeStructures/" + structure1Name;
# StructurePath = "../Scenes/" + sceneName + "/ProbeStructures/" + "0";
print(Structure0Path);
print(Structure1Path);

err0 = 0;
err1 = 0;
BaseCoeffsFilesPath = "C:/Temp/";
FittedCoeffsFilePath = "";
coeffFiles0 = open(BaseCoeffsFilesPath + "guess0.txt", "w");
coeffFiles1 = open(BaseCoeffsFilesPath + "guess1.txt", "w");
coeffFiles2 = open(BaseCoeffsFilesPath + "guess2.txt", "w");

for probeID in xrange(NumberOfProbes):

	coeffFile0 = open(Structure0Path + "/coefficients/Coefficients_" + repr(probeID)+ ".txt", 'r');
	coeffFile1 = open(Structure1Path + "/coefficients/Coefficients_" + repr(probeID)+ ".txt", 'r');
	
	for a in xrange(NumberOfCoefficientsPerProbe):
		r0 = float(coeffFile0.readline());
		g0 = float(coeffFile0.readline());
		b0 = float(coeffFile0.readline());
		r1 = float(coeffFile1.readline());
		g1 = float(coeffFile1.readline());
		b1 = float(coeffFile1.readline());

		err0 += (r0 - r1) ** 2 + (g0 - g1) ** 2 + (b0 - b1) ** 2;

print("err0r: " + repr(err0));