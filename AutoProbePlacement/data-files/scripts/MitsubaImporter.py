#!/usr/bin/python
# Filename: MitsubaImporter.py

import os, sys

sys.path.append('C:/git/mitsuba2/dist/python/3.4/')
# 	# Ensure that Python will be able to find the Mitsuba core libraries
os.environ['PATH'] = 'C:/git/mitsuba2/dist/' + os.pathsep + os.environ['PATH']

print(os.environ['PATH'])
from mitsuba.core import *
from mitsuba.render import SceneHandler, RenderQueue, RenderJob, Scene

# from string import maketrans
import helper