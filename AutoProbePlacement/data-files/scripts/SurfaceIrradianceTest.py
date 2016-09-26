#import mitsuba
import multiprocessing
#import sys
#import os.path

import os, sys
# NOTE: remember to specify paths using FORWARD slashes (i.e. '/' instead of
# '\' to avoid pitfalls with string escaping)
# Configure the search path for the Python extension module
sys.path.append('C:/Users/polardpj.artichaut/Downloads/Mitsuba 0.5.0 32bit/Mitsuba 0.5.0/python/2.7')
# Ensure that Python will be able to find the Mitsuba core libraries
os.environ['PATH'] = 'C:/Users/polardpj.artichaut/Downloads/Mitsuba 0.5.0 32bit/Mitsuba 0.5.0/' + os.pathsep + os.environ['PATH']
import mitsuba

#from mitsuba.core import *
#from mitsuba.render import SceneHandler, RenderQueue, RenderJob

# from string import maketrans


def makeImage(x,y,z):
	scheduler = Scheduler.getInstance()

	#Assign workers to available threads
	for i in range(0, multiprocessing.cpu_count()):
		scheduler.registerWorker(LocalWorker( i, 'wrk%i' % i))
	scheduler.start()
	queue = RenderQueue();
	pmgr = PluginManager.getInstance();
	fileResolver = Thread.getThread().getFileResolver()

	# Link to scene folder
	# fileResolver.appendPath('/Users/joelpp/Documents/Maitrise/Scenes/RedGreenScene/')
	# fileResolver.appendPath('/Users/joelpp/Documents/Maitrise/Scenes/BetterScene/')
	# fileResolver.appendPath('/Users/joelpp/Documents/Maitrise/Scenes/ZScene/')
	fileResolver.appendPath('/Users/joelpp/Documents/Maitrise/Scenes/Z2/')
	# fileResolver.appendPath('/Users/joelpp/Documents/Maitrise/Scenes/BoxTest/')
	# fileResolver.appendPath('/Users/joelpp/Documents/Maitrise/Scenes/BoxPillar/')

	paramMap = StringMap()

	#Start scheduler

	# print(sys.argv[0]);


	queue = RenderQueue()
	total = 9;

	# scene = SceneHandler.loadScene(fileResolver.resolve("RedGreen.xml"), paramMap)
	# scene = SceneHandler.loadScene(fileResolver.resolve("betterScene.xml"), paramMap)
	scene = SceneHandler.loadScene(fileResolver.resolve("Z2.Scene.00051.xml"), paramMap)
	# scene = SceneHandler.loadScene(fileResolver.resolve("box.Scene.00001.xml"), paramMap)


	# Integrator properties
	# integratorType = 'ao'
	integratorType = 'path'
	# integratorType = 'direct'

	# Sampler properties
	sampleCount = 128

	#Create integrator
	integrator = pmgr.create({
		'type' : integratorType,
		})

	#Create sampler
	sampler = pmgr.create({
		'type' : 'ldsampler',
		'sampleCount' : sampleCount
		})

	#Create film
	filmProps = Properties('mfilm')
	filmProps['width'] = 512
	filmProps['height'] = 256
	filmProps['banner'] = False
	filmProps['pixelFormat'] = "rgb"
	film = pmgr.createObject(filmProps)

	offset = 4;
	#Create sensor and transform for sensor
	toWorld = Transform.lookAt(

	    # Point(i-offset, 1+k, j-offset),
	    Point(x, y, z),
	    Point(x+1, y, z),

	    # Point(i-offset, 0, j-offset),
	    # Vector(1, 0, 0)) #up direction before I had to fix stuff
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

	irrSensorProps = Properties('irradiancemeter');
	irrSensor = pmgr.createObject(irrSensorProps);
	irrSensor.addChild('film', film)
	irrSensor.addChild('sampler', sampler)
	irrSensor.configure();
	# scene.getShapes()[4].addChild('sensor', irrSensor);
	# Now, the sensor can be configured
	sensor.configure();
	# sampler.configure();
	scene.setIntegrator(integrator)
	# scene.setSampler(sampler);
	# scene.setSensor(sensor);
	# scene.setSensor(irrSensor);

	scene.configure()
	print(scene.getShapes());

	# Your custom destination goes here.
	destination = '../../Images/IrradianceMeter/'
	destination += 'img.png'
	scene.setDestinationFile(destination)

	#Render!
	job = RenderJob('myRenderJob', scene, queue)
	job.start()
	print(Statistics.getInstance().getStats())

	queue.waitLeft(0)
	queue.join()

#makeImage(0,2,0);

#sys.stdout.write('\a')
#sys.stdout.flush()
# exit();