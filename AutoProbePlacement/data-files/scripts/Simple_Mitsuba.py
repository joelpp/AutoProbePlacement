import multiprocessing
import os, sys
# NOTE: remember to specify paths using FORWARD slashes (i.e. '/' instead of
# '\' to avoid pitfalls with string escaping)
# Configure the search path for the Python extension module
sys.path.append('C:/Users/polardpj.artichaut/Downloads/Mitsuba 0.5.0 32bit/Mitsuba 0.5.0/python/2.7')
# Ensure that Python will be able to find the Mitsuba core libraries
os.environ['PATH'] = 'C:/Users/polardpj.artichaut/Downloads/Mitsuba 0.5.0 32bit/Mitsuba 0.5.0/' + os.pathsep + os.environ['PATH']
import mitsuba

from mitsuba.core import *
from mitsuba.render import RenderQueue, RenderJob
from mitsuba.render import SceneHandler

scheduler = Scheduler.getInstance()
# Start up the scheduling system with one worker per local core
for i in range(0, multiprocessing.cpu_count()):
	scheduler.registerWorker(LocalWorker(i, 'wrk%i' % i))
scheduler.start()



# Get a reference to the thread's file resolver
fileResolver = Thread.getThread().getFileResolver()
# Register any searchs path needed to load scene resources (optional)
fileResolver.appendPath('../ZScene/')
# Optional: supply parameters that can be accessed
# by the scene (e.g. as $myParameter)
paramMap = StringMap()
# Load the scene from an XML file
scene = SceneHandler.loadScene(fileResolver.resolve("Z.Scene.00001.xml"), paramMap)
# Display a textual summary of the scene's contents
print(scene)





# Create a queue for tracking render jobs
queue = RenderQueue()
scene.setDestinationFile('renderedResult')
# Create a render job and insert it into the queue
job = RenderJob('myRenderJob', scene, queue)
job.start()

# Wait for all jobs to finish and release resources
queue.waitLeft(0)
queue.join()
# Print some statistics about the rendering process
print(Statistics.getInstance().getStats())

