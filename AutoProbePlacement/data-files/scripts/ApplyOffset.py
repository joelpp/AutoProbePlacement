import os
import helper
import sys;
(sceneName, structureName) = helper.getSceneAndStructureFromCommandLineArguments(sys.argv);

offsetFile = open("C:/temp/dp.txt", 'r');
logFile = open("C:/temp/dplog.txt", 'a');
probePositionFile = open("../Scenes/" + sceneName + "/ProbeStructures/" + structureName + "/probeList.txt", 'r');
tempPositions = [];

maxDisplacement = 0.25;
counter = 0;
for line in probePositionFile:

	splitLine = line.split();
	basePos = [float(x) for x in splitLine];
	
	dp = [0,0,0];
	dp[0] = float(offsetFile.readline());
	dp[1] = float(offsetFile.readline());
	dp[2] = float(offsetFile.readline());
	# dp[0] *= -1;
	print(dp);
	logFile.write(repr(counter) + " : " + repr(dp) + "\n");
	counter += 1;
	# if ((abs(dp[0]) > abs(dp[1])) and (abs(dp[0]) > abs(dp[2]))):
		# biggest = 0;
		# dp[1] = 0;
		# dp[2] = 0;
	# elif (abs(dp[1]) > abs(dp[2])):
		# biggest = 1;
		# dp[0] = 0;
		# dp[2] = 0;
	# else:
		# biggest = 2;
		# dp[0] = 0;
		# dp[1] = 0;
	
	length = helper.vectorNorm(dp);
	if (length > maxDisplacement):
		print(length);
		factor = length / maxDisplacement;
		print(factor);
		dp = [x / factor for x in dp];
	print(dp);

	newPos = [a-b for a,b in zip(basePos, dp)];
	tempPositions.append(newPos);
probePositionFile.close();

infoFile = open("../Scenes/" + sceneName + "/ProbeStructures/" + structureName + "/info.txt", 'r+');
infoFileContents = infoFile.readlines();
# infoFileContents[0] = "firstProbePosition " + repr(tempPositions[0][0]) + " " + repr(tempPositions[0][1]) + " " + repr(tempPositions[0][2]) + "\n";
infoFile.seek(0);
for line in infoFileContents:
	infoFile.write(line);

probePositionFile = open("../Scenes/" + sceneName + "/ProbeStructures/" + structureName + "/probeList.txt", 'w');
for i in xrange(len(tempPositions)):
	probePositionFile.write(repr(tempPositions[i][0]) + " " + repr(tempPositions[i][1]) + " " + repr(tempPositions[i][2]) + "\n");
probePositionFile.close();

logFile.write("\n\n");
logFile.close();
