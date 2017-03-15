#!/usr/bin/python
# Filename: helper.py

#------------------------------
#|		USEFUL FUNCTIONS	  |	  	
#------------------------------
import math
# import cPickle
# import numpy
from SphericalHarmonics import *
# from sympy.mpmath import *

# from scipy import special
# from scipy import misc


# shArray = cPickle.load(open('/Users/joelpp/Documents/Maitrise/Code/IrradianceMapping/SHvalues/values.p', 'rb'));
print("Loaded SH Array!");
# shImageArray = numpy.zeros((120));
shArray = [];
for i in range(120):
	shArray.append(None);
	
# for i in xrange(120):
# 	shImageArray.append(misc.imread("/Users/joelpp/Documents/Maitrise/Images/IrradianceMaps/SHRender"+repr(i)+".png"));
pi = 3.141592654;
print("Loaded all SH Images!");
	
def getSmallestCoordForInterpolation(position):
	step = 0.5;
	x = math.floor((position[0] / step)) * step;
	y = math.floor((position[1] / step)) * step;
	z = math.floor((position[2] / step)) * step;
	
	return [x,y,z]
def getImgSize(img):
	return (len(img), len(img[0]));

def sphericalToCartesian(r, phi, theta):
	x = r * math.sin(theta) * math.cos(phi);
	y = r * math.sin(theta) * math.sin(phi);
	z = r * math.cos(theta);

	return (x,y,z);

def cartesianToSpherical(x,y,z):
	phi = math.atan2(y, x);

	theta = math.acos(z);

	return (phi, theta);

def sphericalToCartesian_yup(r, phi, theta):
	x = r * math.sin(theta) * math.cos(phi);
	y = r * math.cos(theta);
	z = r * math.sin(theta) * math.sin(phi);

	return (x,y,z);

def cartesianToSpherical_yup(x,y,z):
	phi = math.atan2(z, x);

	theta = math.acos(y);

	return (phi, theta);

#Normalize a vector's components
def normalize(x,y,z):
	length = math.sqrt(x*x +y*y + z*z);
	return (x/length, y/length, z/length);

def normalizeVector(vector):
	# print("Normalizing ["+repr(vector[0])+","+repr(vector[1])+","+repr(vector[2])+",]")
	length = math.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
	return [vector[0]/length, vector[1]/length, vector[2]/length];

def length(x,y,z):
	return math.sqrt(x*x + y*y + z*z);

def vectorNorm(vector):
	return length(vector[0], vector[1], vector[2]);
	
def sinc(x):
	if (math.fabs(x) < 0.0001):
		return 1.0;
	else:
		return (math.sin(x) / x);

def sin(x):
	return math.sin(x);

def show(name, value):
	print(name+" : "+repr(value));

def normalize(x,y,z):
	length = math.sqrt(x*x +y*y + z*z);
	return (x/length, y/length, z/length);

#Transform from a texel location to a normalized U,V coordinate
def ijToUV(i,j, height, width):
	return (j * 2.0 / width, i * 1.0 / height);
	# return ((j-width/2.0)/(width/2.0), (width/2.0 - i)/(width/2.0));

def uvToIJ(u,v, height, width):
	return (height * v, width * u / 2);

#Take in (u,v) coordinates and transform them to (phi, theta) corresponding to an equirectangular projection
def equiUVtoPT(u,v):
	return (math.pi * u, math.pi * v);

def equiPTtoUV(phi, theta):
	pi = 3.14159;
	# print(phi, theta);
	return (phi/pi, theta/pi);

def equiUVtoPT_yup(u,v):

	return (pi*u,pi*(v+1)/2.0); #phi, theta

def ijToUV_yup(i,j,height, width):
	# print(i,j,2.0*j / width -1,2* i / height -1);
	return (2.0*j / width -1, 2.0 * i / height -1);

def uvToIJ_yup(u,v,width,height):
	return (height*(v+1.0)/2.0, width*(u+1.0)/2.0);

def equiPTtoUV_yup(phi, theta):
	return (phi/pi, 2.0*theta/pi - 1.0);

def ijToCartesian_yup(i,j, height, width):
	(u,v) = ijToUV_yup(i,j,height, width);
	(phi, theta) =  equiUVtoPT_yup(u,v);
	return sphericalToCartesian_yup(1, phi, theta);

def cartesianToIJ_yup(x,y,z, height, width):
	(phi, theta) = cartesianToSpherical_yup(x,y,z);
	(u,v) = equiPTtoUV_yup(phi,theta);
	(i,j) = uvToIJ_yup(u,v,width, height);
	return (round(i),round(j));

def phongCoeffs(l, r):
	if (l == 0):
		return pi;
	elif (l == 1):
		return pi * (1.0 + r) / (2.0 + r);
	elif (l == 2):
		return pi * r / (3.0 + r);
	return 0;
	
def NDotOmegaCoeff(l):
	# print(("helper.nDotOmegacvoef",l));
	if (l == 0):
		return math.pi;
	elif (l == 1):
		return 2 * math.pi / 3;
	elif (l == 2):
		return 0.785398;
		# elif(l % 2 == 0):
			# return (2 * math.pi * (-1)**((float(l)/2)-1) * math.factorial(l)) / ((l+2) * (l+1) * 2**l * math.factorial(float(l)/2)**2);
		# else:
			# return 0;
# End of helper.py

def buildStructurePath(prefix, step, dimensions, firstProbePosition):
	return  prefix + "_step_" + repr(step) + \
			"_dimensions_" + repr(dimensions[0]) + "_" + repr(dimensions[1]) + "_" + repr(dimensions[2]) + \
			"_firstProbePosition_" + repr(firstProbePosition[0]) + "_" + repr(firstProbePosition[1]) + "_" + repr(firstProbePosition[2]);

			
def positionFromIndex(index, dimensions, firstProbePosition):
	relativePosition = [index // (dimensions[1] * dimensions[2]), (index // dimensions[2]) % dimensions[1], index % dimensions[2]];
	return [a + b for a, b in zip(firstProbePosition, relativePosition)];
	
def dot(v0, v1):
	return sum([a * b for a, b in zip(v0, v1)]);
	
def intToChannelString(value):
	if (value == 0):
		return "R";
	elif (value == 1):
		return "G";
	elif (value == 2):
		return "B";
	else:
		return "INVALID CHANNEL - IMPLEMENTED VALUES ARE 0,1,2";
		
def buildScenePath(sceneName):
	return "../Scenes/" + sceneName;
		
def buildProbeStructurePath(sceneName, structureName):
	return "../Scenes/" + sceneName + "/ProbeStructures/" + structureName;
	
def getSceneAndStructureFromCommandLineArguments(args):
	sceneName = args[1];
	structureName = args[2];
	
	return (sceneName, structureName);


def getParamFromLine(line):
	splitLine = line.split();
	length = len(splitLine);

	print(splitLine)
	toReturn = [];
	toReturn.append(splitLine[0]);

	splitLine.pop(0);
	
	length = len(splitLine);
	
	if (length == 1):
		toReturn.append(splitLine[0]);
	else:
		toReturn.append([x for x in splitLine]);
		
	return toReturn;
def readSceneInfo(sceneName):
	# define the paths to the info files
	scenePath = "../Scenes/" + sceneName;
	sceneInfoFile = open(scenePath + "/SceneInfo.txt");

	# Put all the lines in lists
	sceneInfo = sceneInfoFile.readlines();

	returnInfo = {};
	seekedInfo = ["scale"];

	for line in sceneInfo:
		for info in seekedInfo:
			if line.startswith(info):
				returnInfo[info] = getParamFromLine(line);

	return returnInfo;
	
	
def readSceneAndProbeStructureInfo(sceneName, probeStructureName):
	
	# define the paths to the info files
	scenePath = "../Scenes/" + sceneName;
	probeStructurePath = scenePath + "/ProbeStructures/" + probeStructureName;
	
	sceneInfoFile = open(scenePath + "/SceneInfo.txt");
	probeStructureInfoFile = open(probeStructurePath + "/info.txt");
	
	# Put all the lines in lists
	sceneInfo = sceneInfoFile.readlines();
	probeStructureInfo = probeStructureInfoFile.readlines();
	
	returnInfo = {};

	# Look for the probe structure type to see what info we'll need
	# for line in probeStructureInfo:
	# 	if line.startswith("type"):
	# 		type = getParamFromLine(line);
	# 		returnInfo[type] = type;

	# if (returnInfo["type"] == "tetrahedral"):
	# 	seekedInfo = ["scale", "dimensions"];
	# else:
	# 	seekedInfo = ["scale", "step", "dimensions", "firstProbePosition", "minBound", "maxBound"];

	# seekedInfo.append("gamma");
	# seekedInfo.append("sampleCount");
	# seekedInfo.append("width");
	# seekedInfo.append("height");
	# seekedInfo.append("integrator");
		
	for line in sceneInfo:
		if (line == "\n"):
			continue;
		param = getParamFromLine(line);
		returnInfo[param[0]] = param[1] 
			
	for line in probeStructureInfo:
		if (line == "\n"):
			continue;
		param = getParamFromLine(line);
		returnInfo[param[0]] = param[1] 

	return returnInfo;

def sRGBToRGBVal(val):
	if (val < 0.04045):
		return val / 12.92;
	else:
		return ((val + 0.055) / 1.055)**2.4; 
	
def sRGBToRGB(srgbColor):
	rgb = [];
	for x in srgbColor:
		rgb.push(sRGBToRGBVal(x))

	return rgb;

def writeVector(file, vector): 
	toWrite = "";
	for i in range(len(vector)):
		toWrite += str(vector[i]);
		toWrite += " ";
	toWrite += "\n";
	file.write(toWrite);

def add(v1, v2):
	return [a+b for a,b in zip(v1, v2)];

def intVector(v):
	return [int(x) for x in v];