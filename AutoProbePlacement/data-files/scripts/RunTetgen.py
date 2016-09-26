import sys
import os
from subprocess import call

if (len(sys.argv) != 3):
	print("Usage: python " + sys.argv[0] + " sceneName probeStructureName");
	exit(1);
sceneName = sys.argv[1];
structureName = sys.argv[2];

if not (os.path.exists( "../Scenes/" + sceneName)):
	print("Scene does not exist!");
	exit(1);

	
call(["C:/libraries/g3d/samples/aTest/data-files/tetgen/Debug/tetgen.exe", "C:/libraries/g3d/samples/aTest/data-files/Scenes/" + sceneName + "/ProbeStructures/" + structureName + "/tetgenOutput/probeList", "-n", "-f", "-e"]);
