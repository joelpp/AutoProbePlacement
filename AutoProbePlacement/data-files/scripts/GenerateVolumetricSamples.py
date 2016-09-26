import random
import onecamera
import helper
from PIL import Image

path = "../temp/volumetricSamples/"
savePath = path+"volumetricSamples.txt";
endl = "\n"
def randomPos(bounds):
	return [random.uniform(bounds['x'][0], bounds['x'][1]), \
			random.uniform(bounds['y'][0], bounds['y'][1]), \
			random.uniform(bounds['z'][0], bounds['z'][1])];

def randomNormal():
	global normalBounds
	return randomPos(normalBounds);


def renderVolumetricSamples():
	pos = randomPos(bounds);
	print pos;

	onecamera.makeProbe(pos[0], pos[1], pos[2], 99,0);

	img = Image.open(path+"99.png");
	print(img)
	normal = randomNormal();
	
	(i,j) = helper.cartesianToIJ_yup(normal[0], normal[1], normal[2], 512, 256);
	rgb = img.getpixel((i,j));
	saveFile = open(savePath, 'a');
	saveFile.write("0"+endl);
	saveFile.write(repr(pos[0]) + " " + repr(pos[1]) + " " + repr(pos[2])+endl);
	saveFile.write("0 0"+endl);
	saveFile.write(repr(normal[0]) + " " + repr(normal[1]) + " " + repr(normal[2])+endl);
	saveFile.write("0 0 0 0"+endl);
	saveFile.write("0 0 0 0"+endl);
	saveFile.write(repr(rgb[0] / 255.0) + " " + repr(rgb[1] / 255.0) + " " + repr(rgb[2] / 255.0)+endl);
	saveFile.write(endl);
	print "(i,j): "+repr((i,j)); 
	print "img[i][j]: "+repr(img.getpixel((i,j)));
	
bounds = {
	'x': [-5, 5],
	'y': [0, 5],
	'z': [-5, 5]
}

normalBounds = {
	'x': [-1, 1],
	'y': [-1, 1],
	'z': [-1, 1]
}

for i in xrange(10):
	renderVolumetricSamples();
