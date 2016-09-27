import os;
import sys;

import helper;

if (len(sys.argv) != 3):
	print("Usage: python " + sys.argv[0] + " sceneName probeStructureName");
	exit(1);
(sceneName, structureName) = helper.getSceneAndStructureFromCommandLineArguments(sys.argv);
globalInfo = helper.readSceneAndProbeStructureInfo(sceneName, structureName);

if not (os.path.exists( "../Scenes/" + sceneName)):
	print("Scene does not exist!");
	exit(1);

rootPath = "../Scenes/" + sceneName + "/" + "ProbeStructures/" + structureName;

probeListFile = open(rootPath + "/probeList.txt", "w+");

if (globalInfo["type"] == "trilinear"):
	initialPos = globalInfo["firstProbePosition"]
	step = float(globalInfo["step"]);
	dim = helper.intVector(globalInfo["dimensions"]);
	print(initialPos, step, dim)
	totalNumber = dim[0] * dim[1] * dim[2];

	for i in xrange(dim[0]):
		for j in xrange(dim[1]):
			for k in xrange(dim[2]):
				pos = helper.add(initialPos, [i*step, j*step, k*step])
				helper.writeVector(probeListFile, pos);