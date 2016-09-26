from PIL import Image

def fixImage(id, commandList):
	pathEnd = repr(id)+".png";
	
	im = Image.open("../temp/surfaceRenders/Reconstruction_optimized/surfaceRender_"+pathEnd);
	for i in xrange(len(commandList)):
		im = im.transpose(commandList[i]);
	im.save("../temp/surfaceRenders/Reconstruction_optimized/surfaceRender_"+pathEnd);
	
	im = Image.open("../temp/surfaceRenders/X/surfaceRenderX_"+pathEnd);
	for i in xrange(len(commandList)):
		im = im.transpose(commandList[i]);
	im.save("../temp/surfaceRenders/X/surfaceRenderX_"+pathEnd);
	
	im = Image.open("../temp/surfaceRenders/Y/surfaceRenderY_"+pathEnd);
	for i in xrange(len(commandList)):
		im = im.transpose(commandList[i]);
	im.save("../temp/surfaceRenders/Y/surfaceRenderY_"+pathEnd);
	
	im = Image.open("../temp/surfaceRenders/Z/surfaceRenderZ_"+pathEnd);
	for i in xrange(len(commandList)):
		im = im.transpose(commandList[i]);
	im.save("../temp/surfaceRenders/Z/surfaceRenderZ_"+pathEnd);


# fixImage(2, [Image.FLIP_LEFT_RIGHT]);
# fixImage(3, [Image.FLIP_TOP_BOTTOM]);
# fixImage(4, [Image.ROTATE_90]);
# fixImage(6, [Image.FLIP_LEFT_RIGHT]);

fixImage(0, [Image.FLIP_TOP_BOTTOM]);
fixImage(1, [Image.FLIP_TOP_BOTTOM]);
fixImage(2, [Image.FLIP_TOP_BOTTOM]);
fixImage(2, [Image.FLIP_LEFT_RIGHT]);
# fixImage(3, [Image.FLIP_TOP_BOTTOM]);
fixImage(4, [Image.FLIP_TOP_BOTTOM]);
fixImage(4, [Image.ROTATE_270]);
fixImage(5, [Image.FLIP_TOP_BOTTOM]);
fixImage(6, [Image.FLIP_TOP_BOTTOM]);
fixImage(7, [Image.FLIP_TOP_BOTTOM]);
fixImage(8, [Image.FLIP_TOP_BOTTOM]);
fixImage(9, [Image.FLIP_TOP_BOTTOM]);
