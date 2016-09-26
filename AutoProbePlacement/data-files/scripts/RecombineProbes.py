import random;
import os;
import helper;


NumberOfCoefficientsPerProbe = 9;
CurrentInfoFile = open("C:/temp/current.txt");
NumberOfReadSamples = int( CurrentInfoFile.readline() );

step = 1;

# sceneName = "ZScene"
# firstProbePosition = [-5,0,-5];
# dimensions = [11,6,11];

# sceneName = "crytek-sponza"
# firstProbePosition = [-10,-1,-7];
# dimensions = [21,10,14];
# firstProbePosition = [1, 4, 0];
# dimensions = [1,1,1];

sceneName = "zcbox"
firstProbePosition = [0,0,0];
dimensions = [7,7,7];

NumberOfProbes = dimensions[0] * dimensions[1] * dimensions[2];
NumberOfElements = NumberOfProbes * NumberOfCoefficientsPerProbe;
print((NumberOfProbes, NumberOfCoefficientsPerProbe))


StructurePath = "../Scenes/" + sceneName + "/ProbeStructures/" + "Fitted_" + repr(NumberOfReadSamples) + "_Samples0";
# StructurePath = "../Scenes/crytek-sponza/ProbeStructures/" + "rgb";

counter = 1;
while os.path.exists(StructurePath):
	StructurePath = StructurePath[:-1] + repr(counter);
	counter += 1;
print(StructurePath);

if not os.path.exists(StructurePath):
    os.makedirs(StructurePath);
    os.makedirs(StructurePath+"/coefficients/");
    os.makedirs(StructurePath+"/projhbes/");
	
infoFile = open(StructurePath + "/info.txt", "w");
infoFile.write("step " + repr(step) + "\n");
infoFile.write("dimensions " + repr(dimensions[0]) + " " +  repr(dimensions[1]) + " " + repr(dimensions[2]) + "\n");
infoFile.write("firstProbePosition " + repr(firstProbePosition[0]) + " " +  repr(firstProbePosition[1]) + " " + repr(firstProbePosition[2]) + "\n");

probeListFile = open(StructurePath + "/probeList.txt", "w");

BaseCoeffsFilesPath = "C:/Temp/";
FittedCoeffsFilePath = "";
coeffFiles0 = open(BaseCoeffsFilesPath + "v0.txt", "r");
coeffFiles1 = open(BaseCoeffsFilesPath + "v1.txt", "r");
coeffFiles2 = open(BaseCoeffsFilesPath + "v2.txt", "r");

for probeID in xrange(NumberOfProbes):
	
	probePosition = helper.positionFromIndex(probeID, dimensions, firstProbePosition);
	
	probeListFile.write(repr(probePosition[0]) + " " + repr(probePosition[1]) + " " + repr(probePosition[2]) + "\n");
	
	coeffFile = open(StructurePath + "/coefficients/Coefficients_" + repr(probeID)+ ".txt", 'w');
	
	for i in xrange(NumberOfCoefficientsPerProbe):
		coeffFile.write(repr( float(coeffFiles0.readline() ) ) + "\n")
		coeffFile.write(repr( float(coeffFiles1.readline() ) ) + "\n")
		coeffFile.write(repr( float(coeffFiles2.readline() ) ) + "\n")
		
