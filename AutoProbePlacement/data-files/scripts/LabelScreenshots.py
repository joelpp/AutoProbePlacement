import sys
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw

def label(sceneName, ID):
	root = "../Scenes/" + sceneName + "/Optimizations/" + ID + "/";
	errorLogFile = open(root + "errorlog.txt",'r');

	errorlines = errorLogFile.readlines();

	for i in xrange(len(errorlines)):
		
		try:
			path = root + "screens/" + repr(i+1) + ".jpg"
			
			image = Image.open(path);
			draw = ImageDraw.Draw(image)
			
			
			error = repr(i) + " " + errorlines[i];
			
			font = ImageFont.truetype("arial.ttf", 16)
			draw.text((0, 0), error, (255,255,255), font=font)
			image.save(path)
		except:
			pass;
			
	errorLogFile.close();