from PIL import Image;
import numpy as np;
import sys;
# print(len(new_img));
# for i in xrange(int(sys.argv[1])):
	# renderResult = Image.open('../crytek-sponza/irradiance/irr_'+repr(i)+'.png');
	# print((i, list(renderResult.getdata())));
renderResult = Image.open('../crytek-sponza/MASKTEST.png');
print(list(renderResult.getdata()));