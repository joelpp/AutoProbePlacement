import os;
import sys;

sceneName = sys.argv[1];
structureNames = [sys.argv[2], sys.argv[3]];
probeIDs = ["0", "0"];

pathToProbes = ["../Scenes/" + sceneName + "/ProbeStructures/" + x + "/Coefficients/Coefficients_" for x in structureNames];

probePaths = [y + x + ".txt" for x,y in zip(probeIDs, pathToProbes)];

probeFiles = [open(x, 'r') for x in probePaths]

probeFileContents = [x.readlines() for x in probeFiles];

difference = [float(x) - float(y) for x,y in zip(probeFileContents[0], probeFileContents[1])]


c = 0;
i = 0;
for line in difference:
	print([i, c, line])

	c += 1;
	if (c == 3):
		c = 0;
		i += 1;

