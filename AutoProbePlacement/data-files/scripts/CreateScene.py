import os;
import sys;

if (len(sys.argv) != 2):
	print("Usage: python " + sys.argv[0] + " sceneName");
	exit(1);
sceneName = sys.argv[1];
rootPath = "../Scenes/" + sceneName;
if not (os.path.exists(rootPath)):
	
	os.mkdir(rootPath);
	os.mkdir(rootPath + "/ProbeStructures");
	os.mkdir(rootPath + "/SampleSets");
	os.mkdir(rootPath + "/objs");
	open(rootPath + "MitsubaScene.xml", 'w')
	open(rootPath + "SceneInfo.txt", 'w')