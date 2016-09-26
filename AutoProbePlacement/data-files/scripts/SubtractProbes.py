import os;
import sys;

sceneName = sys.argv[1];
structureName = sys.argv[2];
probeIDs = [sys.argv[3], sys.argv[4]];

pathToProbes = "../Scenes/" + sceneName + "/ProbeStructures/" + structureName + "/Coefficients/Coefficients_";

probePaths = [pathToProbes + x + ".txt" for x in probeIDs];

probeFiles = [open(x, 'r') for x in probePaths]

probeFileContents = [x.readlines() for x in probeFiles];
print(probeFileContents);
print(float(probeFileContents[0][0]));
difference = [float(x) - float(y) for x,y in zip(probeFileContents[0], probeFileContents[1])]

for line in difference:
	print(line)