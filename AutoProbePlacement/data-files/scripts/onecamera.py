import multiprocessing
import os, sys
# NOTE: remember to specify paths using FORWARD slashes (i.e. '/' instead of
# '\' to avoid pitfalls with string escaping)
# # Configure the search path for the Python extension module
# # Ensure that Python will be able to find the Mitsuba core libraries

# Ensure that Python will be able to find the Mitsuba core libraries

if (True):
	sys.path.append('C:/git/mitsuba/dist/python/2.7/')
	# Ensure that Python will be able to find the Mitsuba core libraries
	os.environ['PATH'] = 'C:/git/mitsuba/dist/' + os.pathsep + os.environ['PATH']
else:
	sys.path.append('C:/Users/Joel/Downloads/mitsuba-c7aac473729a/mitsuba-c7aac473729a/dist/python/2.7/')
	os.environ['PATH'] = 'C:/Users/Joel/Downloads/mitsuba-c7aac473729a/mitsuba-c7aac473729a/dist/' + os.pathsep + os.environ['PATH']


from mitsuba.core import *
from mitsuba.render import SceneHandler, RenderQueue, RenderJob

from string import maketrans
import helper

# def readProbeStructureInfo():
(sceneName, structureName) = helper.getSceneAndStructureFromCommandLineArguments(sys.argv);
globalInfo = helper.readSceneAndProbeStructureInfo(sceneName, structureName);

pmgr = PluginManager.getInstance();
fileResolver = Thread.getThread().getFileResolver();

paramMap = StringMap()
path = "C:/git/AutoProbePlacement/AutoProbePlacement/Scenes/" + sceneName;
fileResolver.appendPath(path);


renderType = sys.argv[3];


scene = SceneHandler.loadScene("../Scenes/zcbox/MitsubaScene2.xml", paramMap)

#Initialise Mitsuba stuff
queue = RenderQueue()

scheduler = Scheduler.getInstance()

#Assign workers to available threads
for i in range(0, multiprocessing.cpu_count()):
	scheduler.registerWorker(LocalWorker( i, 'wrk%i' % i))
scheduler.start()

def makeProbe(x, y, z, probeCount, rootPath, pRenderType):

	total = 9;

	# Integrator properties
	integratorType = pRenderType[:-1].lower();

	# Sampler properties
	sampleCount = 1
	if (pRenderType == "Probes"):
		sampleCount = 128;
		integratorType = "path";
		# integratorType = "custom_path";

	#Create integrator
	integrator = pmgr.create({
		'type' : integratorType,
		'maxDist': 900.0,
		'samplePosition': Point(0,0,0),
		'sampleNormal': Normal(0,0,0),
		'pointSampling': 0,
		'logFinalValue': False,
		'logRayTerminations': False
		})

	#Create sampler
	sampler = pmgr.create({
		'type' : 'ldsampler',
		'sampleCount' : sampleCount
		})

	#Create film
	filmProps = Properties('hdrfilm')
	if (pRenderType == "Probes"):
		filmProps = Properties('ldrfilm')
	filmProps['width'] = 128
	filmProps['height'] = 64	
	filmProps['banner'] = False
	filmProps['pixelFormat'] = "rgb"
	# filmProps['gamma'] = 0.0
	filmProps['rfilter'] = "box"
	filmProps['componentFormat'] = "float32"
	
	film = pmgr.createObject(filmProps)
	
	# boxProps = Properties('box');
	# boxProps['radius'] = 0;
	# film.setFilter(pmgr.createObject(boxProps));

	# offset = 4;
	scale = globalInfo["scale"];
	x /= scale;
	y /= scale;
	z /= scale;
	
	#Create sensor and transform for sensor
	toWorld = Transform.lookAt(

	    Point(x, y, z),
	    Point(x+1, y, z),
	    Vector(0, 1, 0)) #up direction after I fixed stuff
		
	sensorProps = Properties('spherical')
	# sensorProps = Properties('perspective')
	sensorProps['toWorld'] = toWorld
	sensorProps['fov'] = 45.0
	sensor = pmgr.createObject(sensorProps)
	#format of a lookat transform: lookAt(origin, target, updirection)

	#Make it all come together
	film.configure() 
	sensor.addChild('film', film)
	sensor.addChild('sampler', sampler)


	# Now, the sensor can be configured
	sensor.configure();
	# sampler.configure();
	scene.setIntegrator(integrator)
	scene.setSampler(sampler);
	scene.addSensor(sensor);
	scene.setSensor(sensor);


	scene.configure()
	# print(scene);
	
	# Your custom destination goes here.
	destination = rootPath + "/" + pRenderType + "/" + pRenderType[:-1] + "_"
	destination += repr(probeCount)+'.png'
	scene.setDestinationFile(destination)
	print(destination)
	#Render!
	job = RenderJob('myRenderJob', scene, queue)
	job.start()
	# print(Statistics.getInstance().getStats())

	queue.waitLeft(0)
	queue.join()
	
def createProbeStructure(rootPath):
	
	if not (os.path.exists(rootPath)):
		os.mkdir(rootPath);
		os.mkdir(rootPath + "/probes");
		os.mkdir(rootPath + "/coefficients");
		os.mkdir(rootPath + "/renders");
		os.mkdir(rootPath + "/misc");
	if not (os.path.exists(rootPath + "/" + renderType)):
		os.mkdir(rootPath + "/" + renderType);

def makeProbeGrid(rootPath, _offset, _step, _numIterations):

	
	fileList = open(rootPath + '/probeList.txt', 'r');

		
	renderList = [];
	probeCount = 0;
	for i in xrange(int(_numIterations)):
		# for j in xrange(_numIterations[1]):
			# for k in xrange(_numIterations[2]):

				# x = i*_step[0]+_offset[0];
				# y = j*_step[1]+_offset[1];
				# z = k*_step[2]+_offset[2];
		splitLine = fileList.readline().split();
		print(splitLine);
		(x, y, z) = (float(x) for x in splitLine);
		print((x,y,z));
		# exit();
		renderList.append([x, y, z]);
		# fileList.write(repr(x) + ' ' + repr(y) + ' ' + repr(z) + '\r\n');
		print("Rendering probe #" + repr(probeCount) + ", pos: " + repr((x,y,z)));
		makeProbe(x,y,z,probeCount, rootPath, renderType);
		probeCount += 1;
				

def makeProbeList(rootPath):
	fileList = open(rootPath + '/probeList.txt', 'r');
	probeCount = 0;
	for line in fileList:
		splitLine = line.split();
		print(splitLine);
		(x, y, z) = (float(x) for x in splitLine);
		print((x,y,z));
		
		if (renderType == "all"):
			makeProbe(x,y,z,probeCount, rootPath, "Probes");
			makeProbe(x,y,z,probeCount, rootPath, "Positions");
			makeProbe(x,y,z,probeCount, rootPath, "Normals");
		else:
			makeProbe(x,y,z,probeCount, rootPath, renderType);
		probeCount += 1;


def makeProbeAt(index, position):
	x = position[0];
	y = position[1];
	z = position[2];
	makeProbe(x,y,z,index, 0, renderType);


def renderProbe(rootPath, id):
	fileList = open(rootPath + '/probeList.txt', 'r');
	
	lines = fileList.readlines();
	splitLine = lines[id].split();
	(x, y, z) = (float(x) for x in splitLine);
	makeProbe(x, y, z, id, rootPath, renderType);
	


rootPath = "../Scenes/" + sceneName + "/ProbeStructures/" + structureName;
print(rootPath);
	

createProbeStructure(rootPath);

# If we want to render a single probe
if (len(sys.argv) == 5):
	id = int(sys.argv[4]);
	renderProbe(rootPath, id);
else:
# elif (globalInfo["type"] != "trilinear"):
	makeProbeList(rootPath);
# else:
	# makeProbeGrid(rootPath, globalInfo["firstProbePosition"], globalInfo["step"], globalInfo["dimensions"]);

file = open("C:/temp/log2.txt", 'w');
file.write("EVEN MORE LOL");
file.close();