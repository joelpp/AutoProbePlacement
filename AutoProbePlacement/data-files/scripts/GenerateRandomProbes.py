import random;
import os;
import helper;


NumberOfCoefficientsPerProbe = 9;

step = 1;


# firstProbePosition = [-5,0,-5];
# dimensions = [11,6,11];

# firstProbePosition = [-10,-1,-7];
# dimensions = [21,10,14];

# firstProbePosition = [1, 4, 0];
# dimensions = [1,1,1];

sceneName = "zcbox";
firstProbePosition = [0, 0, 0];
dimensions = [7, 7, 7];

NumberOfProbes = dimensions[0] * dimensions[1] * dimensions[2];
NumberOfElements = NumberOfProbes * NumberOfCoefficientsPerProbe;

StructurePath = "../Scenes/" + sceneName + "/ProbeStructures/" + helper.buildStructurePath("RandomCoefficients", step, dimensions, firstProbePosition);
print(StructurePath);

if not os.path.exists(StructurePath):
    os.makedirs(StructurePath);
    os.makedirs(StructurePath+"/coefficients/");
    os.makedirs(StructurePath+"/probes/");
	
infoFile = open(StructurePath + "/info.txt", "w");
infoFile.write("step " + repr(step) + "\n");
infoFile.write("dimensions " + repr(dimensions[0]) + " " +  repr(dimensions[1]) + " " + repr(dimensions[2]) + "\n");
infoFile.write("firstProbePosition " + repr(firstProbePosition[0]) + " " +  repr(firstProbePosition[1]) + " " + repr(firstProbePosition[2]) + "\n");

probeListFile = open(StructurePath + "/probeList.txt", "w");

for probeID in xrange(NumberOfProbes):
	
	probePosition = helper.positionFromIndex(probeID, dimensions, firstProbePosition);
	
	probeListFile.write(repr(probePosition[0]) + " " + repr(probePosition[1]) + " " + repr(probePosition[2]) + "\n");
	
	coeffFile = open(StructurePath + "/coefficients/Coefficients_" + repr(probeID)+ ".txt", 'w');
	
	for i in xrange(NumberOfCoefficientsPerProbe * 3):
		coeffFile.write(repr(random.uniform(0,1)) + "\n")