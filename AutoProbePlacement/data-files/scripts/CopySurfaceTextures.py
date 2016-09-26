import shutil

path = "../temp/surfaceRenders/"

for i in xrange(9):
	shutil.copy(path+"surfaceRender_"+repr(i)+".png", path+"Reconstruction_optimized/")
	shutil.copy(path+"surfaceRenderReference_"+repr(i)+".png", path+"Reference/")
	shutil.copy(path+"surfaceRenderX_"+repr(i)+".png", path+"X/")
	shutil.copy(path+"surfaceRenderY_"+repr(i)+".png", path+"Y/")
	shutil.copy(path+"surfaceRenderZ_"+repr(i)+".png", path+"Z/")
	