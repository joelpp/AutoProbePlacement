import helper;
import sys;
import scipy.misc as misc;


def markImage(img, _i, _j):
	i = int(_i)-1;
	j = int(_j)-1;
	# print(i,j);
	# print(len(img), len(img[0]))
	img[i][j] = 255;
	# img[i][j][1] = 0;
	# img[i][j][2] = 0;
	return img;

def getIJ():
	a = open("../Validation/UV/uv.txt");

	# v = float(a.readline())*2-1;
	# u = float(a.readline())*2-1;
	# u = 0;
	# v = 0;
	# (i,j) = helper.uvToIJ_yup(u,v, 256, 512);
	# 
	x = float(a.readline());
	y = float(a.readline());
	z = float(a.readline());

	(i,j) = helper.cartesianToIJ_yup(x,y,z, 256, 512);
	# print (u,v);
	print (x,y,z);
	print (i,j);
	return (i,j);

(j,i) = getIJ();
# img = misc.imread("/Users/joelpp/Documents/Maitrise/Scenes/Z2/objparts/largegreenwall_bake.png");
img = misc.imread("/Users/joelpp/Documents/Maitrise/Code/Validation/512x256/SHrender0.png");
img = markImage(img, i, j);
misc.imsave("/Users/joelpp/Documents/Maitrise/Images/IrradianceMeter/normal.png",img);