import sys;
import os;

sceneName = sys.argv[1];
structureName = sys.argv[2];

probeListFile = open("../Scenes/" + sceneName + "/ProbeStructures/" + structureName + "/probeList.txt", 'r');
nodeFile = open("../Scenes/" + sceneName + "/ProbeStructures/" + structureName + "/tetgenOutput/probeList.node", 'w+');

probeListContents = probeListFile.readlines();

nodeFile.write(repr(len(probeListContents)) + " 3 0 0\n");

for i in xrange(len(probeListContents)):
	nodeFile.write(repr(i+1) + " " + probeListContents[i]);