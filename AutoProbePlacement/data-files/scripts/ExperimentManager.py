import shutil
import os
import onecamera
import time


rootPath = "/Users/joelpp/Documents/Maitrise/Code/experiments/";

def setCurrent(i):
	# experimentInfo = open(rootPath + repr(i) + "/info.txt");

	# line = experimentInfo.readline();

	# probeCount = int(line.split(" ")[1]);

	shutil.rmtree(rootPath+"current");
	# shutil.rmtree(rootPath+"current/probes");
	# shutil.rmtree(rootPath+"current/coefficients");

	# os.mkdir(rootPath+"current/probes");
	# os.mkdir(rootPath+"current/coefficients");

	# shutil.copy(rootPath+repr(i)+"/probeList.txt", rootPath+"current");
	# shutil.copy(rootPath+repr(i)+"/info.txt", rootPath+"current");
	shutil.copytree(rootPath+repr(i), rootPath+"current");

def getValueFromFile(_file, _line):

	myFile = open(rootPath+_file);

	for line in myFile:
		splitLine = line.split(" ");
		if splitLine[0] == _line:
			splitLine.pop(0);
			return splitLine;

	return None;

def setValueInFile(_file, _line, _value):
	newFile = open(_file,"w+");

	myFile = open(rootPath+_file, "r+");
	for line in myFile:
		splitLine = line.split(" ");
		if splitLine[0] == _line:
			newFile.write(splitLine[0] + " " + repr(_value)+" \n");
		else:
			newFile.write(line);

	myFile.close();
	newFile.close();
	shutil.copy(_file, rootPath+_file);

	os.remove(_file);


def makeNewExperiment():
	experimentID = int(getValueFromFile("ManagerInfo.txt", "numberOfExperiments")[0]);
	# experimentID = 48
	os.mkdir(rootPath+repr(experimentID));
	os.mkdir(rootPath+repr(experimentID)+"/probes");
	os.mkdir(rootPath+repr(experimentID)+"/coefficients");
	os.mkdir(rootPath+repr(experimentID)+"/renders");
	os.mkdir(rootPath+repr(experimentID)+"/misc");

	probeList = open(rootPath+repr(experimentID)+"/probeList.txt","w+");
	info = open(rootPath+repr(experimentID)+"/info.txt","w+");
	offset = [-5,0,-5];
	step = [1, 1, 1];
	# numIterations = [1,1,1];
	numIterations = [11,5,11];

	numberOfProbes = numIterations[0] * numIterations[1] * numIterations[2];

	onecamera.makeProbeGrid(experimentID, offset, step, numIterations);

	# renderList = [[-1.5,1,4.5],[4.5,1,4.5],[4.5,1,3.5],[1.5,1,3.5],[1.5,1,-4.5],[-1.5,1,-3.5],[-1.5,1,-3.5],[-4.5,1,-3.5],[-4.5,1,-4.5], [-1.5,4,4.5],[4.5,4,4.5],[4.5,4,3.5],[1.5,4,3.5],[1.5,4,-4.5],[-1.5,4,-3.5],[-4.5,4,-3.5],[-4.5,4,-4.5],[0,2.5,0],[0,2.5,4],[0,2.5,-4],[4,2.5,4],[-4,2.5,-4]]; #
	# renderList = [[0, 2.5, 0],[-1.5, 1, -4],[-1.5, 1, 4],[1.5, 1, -4],[1.5, 1, 4],[-2, 4, -4],[-1.5, 4, 4],[1.5, 4, -4],[2, 4, 4]];
	# renderList = [[4.5, 1, 4.5], [4.5,4,4.5], [4.5,1,3.5], [4.5,4,3.5], [1.5,1,3.5], [1.5,4,3.5], [1.5, 1, 0],[1.5,4,0],[-1.5, 1, 4.5],[-1.5,4,4.5],[-1.5, 1, 0],[-1.5,4,0]];
	# renderList = [[0,3,0]];
	# numberOfProbes = len(renderList); onecamera.makeProbeList(experimentID, renderList);

	info.write("numberOfProbes "+repr(numberOfProbes)+" \n");
	# info.write("offset "+repr(offset[0])+" "+repr(offset[1])+" "+repr(offset[2])+"\n");
	# info.write("step "+repr(step[0])+" "+repr(step[1])+" "+repr(step[2])+"\n");

	setValueInFile("ManagerInfo.txt", "numberOfExperiments", experimentID+1);

def addProbes(experimentID, renderList):
	# numberOfProbes = int(getValueFromFile(repr(experimentID)+"/info.txt", "numberOfProbes")[0]);
	# print("number of probes before : "+repr(numberOfProbes));

	onecamera.makeProbeList(experimentID,renderList);
	# setValueInFile(repr(experimentID) + "/info.txt", "numberOfProbes", numberOfProbes+1);


def printExperimentInfo(n):
	print("Experiment "+repr(n)+":");
	path, dirs, files = os.walk(rootPath+repr(n)+"/probes/").next();

	numberOfProbes = len(files);
	print("Number of Probes: "+repr(numberOfProbes));


def printAllExpInfo():
	path, dirs, files = os.walk(rootPath).next();
	numberOfExperiments = len(dirs)-1;
	# print(numberOfExperiments);
	for i in xrange(numberOfExperiments):
		printExperimentInfo(i);
		# print;



# makeNewExperiment()
setCurrent(44);
# printExperimentInfo(10);
# printAllExpInfo();
# addProbes(31, [[-1.5,3,0]]);
