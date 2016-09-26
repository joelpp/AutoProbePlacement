import sys
import helper

def LineToVectorF(line):
	vector = [float(x) for x in line.split()]
	# print(vector);
	return (vector);
	
	
(sceneName, sampleSetName) = helper.getSceneAndStructureFromCommandLineArguments(sys.argv);
NumberOfSamples = int(sys.argv[3]);
globalInfo = helper.readSceneInfo(sceneName);

scenePath = "../Scenes/" + sceneName;
sampleSetPath = scenePath + "/SampleSets/" + sampleSetName; 
samplesFilePath = sampleSetPath + "/IrradianceResults2.txt";
samplesFile = open(samplesFilePath, 'r');
refFile = open('C:/temp/samplesRGB_ref.txt', 'w+');
lines = samplesFile.readlines();


# for linenumber in xrange(len(lines)):
for linenumber in xrange(NumberOfSamples):
	line = lines[linenumber];
	v3 = LineToVectorF(line);
	
	for i in xrange(3):
		refFile.write(repr(v3[i] / 255.0) + "\n");