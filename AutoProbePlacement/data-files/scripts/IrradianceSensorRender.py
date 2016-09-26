import multiprocessing
import os, sys
from PIL import Image
import imageio
import helper

if (True):
	sys.path.append('C:/Users/polardpj.artichaut/Downloads/mitsuba-eaff1cd989f3/mitsuba-eaff1cd989f3/dist/python/2.7/')
	# Ensure that Python will be able to find the Mitsuba core libraries
	os.environ['PATH'] = 'C:/Users/polardpj.artichaut/Downloads/mitsuba-eaff1cd989f3/mitsuba-eaff1cd989f3/dist/' + os.pathsep + os.environ['PATH']
else:
	sys.path.append('C:/Users/Joel/Downloads/mitsuba-c7aac473729a/mitsuba-c7aac473729a/dist/python/2.7/')
	os.environ['PATH'] = 'C:/Users/Joel/Downloads/mitsuba-c7aac473729a/mitsuba-c7aac473729a/dist/' + os.pathsep + os.environ['PATH']

from mitsuba.core import *
from mitsuba.render import SceneHandler, RenderQueue, RenderJob, Scene, Integrator

(sceneName, sampleSetName) = helper.getSceneAndStructureFromCommandLineArguments(sys.argv);
numberOfSamples = int(sys.argv[3]);

globalInfo = helper.readSceneInfo(sceneName);

scenePath = "../Scenes/" + sceneName;
sampleSetPath = scenePath + "/SampleSets/" + sampleSetName; 
samplesFilePath = sampleSetPath + "/SamplePositions.txt";
samplesFile = open(samplesFilePath, 'r');
	
LinesPerSample = 9;

#Initialise Mitsuba stuff
pmgr = PluginManager.getInstance();
fileResolver = Thread.getThread().getFileResolver()

# Link to scene folder
fileResolver.appendPath(scenePath);

paramMap = StringMap()
queue = RenderQueue()
scheduler = Scheduler.getInstance()

sceneScale = globalInfo["scale"];
originalScene = SceneHandler.loadScene(fileResolver.resolve("MitsubaScene.xml"), paramMap)

originalScene.initialize();

# Integrator properties
integratorType = 'custom_path'
# integratorType = 'path'
# integratorType = 'direct'


#Assign workers to available threads
for i in range(0, multiprocessing.cpu_count()):
	scheduler.registerWorker(LocalWorker( i, 'wrk%i' % i))
scheduler.start()

def LineToVectorF(line):
	vector = [float(x) for x in line.split()]
	# print(vector);
	return (vector);

def getNextSample():

	samplesFile.readline(); #scene sample text
	samplesFile.readline(); # #position text

	line = samplesFile.readline(); #position
	position = LineToVectorF(line);

	samplesFile.readline(); # #normal text
	# samplesFile.readline(); #UV

	line = samplesFile.readline(); #Normal
	normal = LineToVectorF(line);
	# samplesFile.readline(); #blank (fossil from tetrahadral days)
	# samplesFile.readline(); #blank (fossil from tetrahadral days)
	# samplesFile.readline(); #diffuse reflectance
	samplesFile.readline(); #empty line
	return (position, normal);

def render(sampleNumber):


	sceneCopy = Scene(originalScene)
	(position, normal) = getNextSample();
	print((position, normal));
	#Create integrator
	integrator = pmgr.create({
		'type' : integratorType,
		'samplePosition': Point(position[0], position[1], position[2]),
		'sampleNormal': Normal(normal[0], normal[1], normal[2]),
		'pointSampling': 1,
		'logFinalValue': False,
		'logRayTerminations': False
		});
	integrator.configure();
	sceneCopy.setIntegrator(integrator);
	
	sampleCount = 8192;
	sampler = pmgr.create({
	'type' : 'ldsampler',
	'sampleCount' : sampleCount
	})

	#Create film
	filmProps = Properties('ldrfilm')
	filmProps['width'] = 1
	filmProps['height'] = 1	
	filmProps['banner'] = False
	filmProps['pixelFormat'] = "rgb"
	# filmProps['gamma'] = 0.0
	filmProps['rfilter'] = "box"
	filmProps['componentFormat'] = "float32"
	film = pmgr.createObject(filmProps)
	film.configure();

	# Now, the sensor can be configured
	sensorProps = Properties('perspective')
	sensorProps['fov'] = 45.0
	sensor = pmgr.createObject(sensorProps)
	
	sensor.addChild('film', film)
	sensor.addChild('sampler', sampler);
	sensor.configure();
	
	sceneCopy.addSensor(sensor);
	sceneCopy.setSensor(sensor);
	
	sceneCopy.setDestinationFile(sampleSetPath + '/irradiance/irr_' + repr(sampleNumber) + '.png');
	
	sceneCopy.configure()
	# print(sceneCopy);

	#Render!
	job = RenderJob('myRenderJob', sceneCopy, queue);
	job.start()
	# print(Statistics.getInstance().getStats())

	queue.waitLeft(0)
	queue.join()

####################
# START OF PROGRAM #
####################


for sampleNumber in xrange(numberOfSamples):
	render(sampleNumber);
	pass;
	
for sampleNumber in xrange(numberOfSamples):
	renderResult = imageio.imread(sampleSetPath + '/irradiance/irr_' + repr(sampleNumber) + '.png');
	resultsFile = open(sampleSetPath + '/IrradianceResults2.txt', 'a');
	result = renderResult[0][0];
	print result;
	resultsFile.write( repr(result[0]) + " " + repr(result[1]) + " " + repr(result[2]) + "\n");

samplesFile.close();