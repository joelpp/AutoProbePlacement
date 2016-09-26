import os;
import sys;
import helper;


(SceneName, SampleSetName) = helper.getSceneAndStructureFromCommandLineArguments(sys.argv);
NumberOfSamples = int(sys.argv[3]);


RGBFile = open("C:/temp/samplesRGB.txt", 'r');

targetFilePath = "../Scenes/" + SceneName + "/SampleSets/" + SampleSetName + "/IrradianceResults2.txt";
targetFile = open(targetFilePath, 'w');

for sampleID in xrange(NumberOfSamples):
	r = float(RGBFile.readline()) * 255;
	g = float(RGBFile.readline()) * 255;
	b = float(RGBFile.readline()) * 255;
	
	targetFile.write(repr(r) + " " + repr(g) + " " + repr(b) + "\n");
	