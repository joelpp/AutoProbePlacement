import os;
import sys;

if (len(sys.argv) != 3):
	print("Usage: python " + sys.argv[0] + " sceneName probeStructureName");
	exit(1);
sceneName = sys.argv[1];
structureName = sys.argv[2];

if not (os.path.exists( "../Scenes/" + sceneName)):
	print("Scene does not exist!");
	exit(1);

rootPath = "../Scenes/" + sceneName + "/" + "ProbeStructures/" + structureName;
if not (os.path.exists(rootPath)):
	os.mkdir(rootPath);
	os.mkdir(rootPath + "/Probes");
	os.mkdir(rootPath + "/Normals");
	os.mkdir(rootPath + "/Positions");
	os.mkdir(rootPath + "/Depths");
	os.mkdir(rootPath + "/Coefficients");
	open(rootPath + "/info.txt", "w+");
	open(rootPath + "/probeList.txt", "w+");