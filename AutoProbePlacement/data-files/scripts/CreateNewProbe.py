
import helper
import sys

if (len(sys.argv) < 8):
	dpFile = open("C:/temp/dp.txt");
	lines = dpFile.readlines();
	displacement = [float(x) for x in lines];
	print(displacement);
	# exit(1);
	
else:
	displacement = [float(sys.argv[5]),float(sys.argv[6]), float(sys.argv[7])];
(sceneName, structureName) = helper.getSceneAndStructureFromCommandLineArguments(sys.argv);
globalInfo = helper.readSceneAndProbeStructureInfo(sceneName, structureName);
sourceProbeID = int(sys.argv[3]);
newProbeID = int(sys.argv[4]);
# displacement = [float(sys.argv[5]),float(sys.argv[6]), float(sys.argv[7])];

probeStructurePath = helper.buildProbeStructurePath(sceneName, structureName);
probeFilePath = probeStructurePath + "/Coefficients/Coefficients_" + repr(sourceProbeID) + ".txt";
probeGradientFilePath = probeStructurePath + "/Coefficients/CoefficientsGradients_" + repr(sourceProbeID) + ".txt";
newProbePath = probeStructurePath + "/Coefficients/Coefficients_" + repr(newProbeID) + ".txt";

probeFile = open(probeFilePath, 'r');
probeGradientFile = open(probeGradientFilePath, 'r');

coeffs = [];
gradients = [];
for line in probeFile:
	coeffs.append(float(line));
probeFile.close();
outputFile = open(newProbePath, 'w+');

axis = 0;
counter = 0;
dot = 0;
for line in probeGradientFile:
	if (line[0] == '/'):
		continue;
	else:
		dot += float(line) * displacement[axis];
		
		axis = (axis + 1);
		if (axis == 3):	
			coeffs[counter] += dot;
			dot = 0;
			axis = 0;
			counter += 1;
		
for coeff in coeffs:
	outputFile.write(repr(coeff) + "\n");
	
outputFile.close();
probeGradientFile.close();