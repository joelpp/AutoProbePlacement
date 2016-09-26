import helper; 
from scipy import misc;
import os;
import sys;
import math;

def getGradient(img):
	
	# (height, width) = helper.getImgSize(img);
	(height, width) = (128,256);

	# gradient = [0,3.14];
	# gradient = [-0.713,-0.535,0.452];
	gradient = [0,0,0];
	luminanceSum = 0;

	for i in range(height):
		for j in range(width):

			(u,v) = helper.ijToUV_yup(i,j, height, width);
			
			(phi, theta) = helper.equiUVtoPT_yup(u,v);

			#Convert from theta,phi to x,y,z coordinates
			(x,y,z) = helper.sphericalToCartesian_yup(1,phi, theta);
			# print(phi, theta);
			R = float(img[i][j][0]);
			G = float(img[i][j][1]);
			B = float(img[i][j][2]);
			
			# if ((R,G,B) != (0,0,0)):
			# 	print("Non-zero pixel at (i,j) = ("+repr(i)+","+repr(j)+"),(u,v) = ("+repr(u)+","+repr(v)+"), (x,y,z) = ("+repr(x)+","+repr(y)+","+repr(z)+")");
			# 	print("Gradient now: "+repr(gradient));
			luminance = 0.2126 * R + 0.7152 * G + 0.0722 * B;
			# luminance =  R + G + B;
			# luminanceSum = luminanceSum + luminance;
			# luminance = 1;
			luminance *= math.sin(theta);
	# 		# gradient[0] += u * luminance;
	# 		# gradient[1] += v * luminance;
			
			gradient[0] += x * luminance;
			gradient[1] += y * luminance;
			gradient[2] += z * luminance;
			# print(gradient);
			
	# # (phi, theta) = helper.equiUVtoPT_yup(gradient[0], gradient[1]);

	# #Convert from theta,phi to x,y,z coordinates
	# # (x,y,z) = helper.sphericalToCartesian_yup(1, phi, theta);
	epsilon = 0.0000001;
	# print("Gradient calculation before normalization or cancellation: "+repr(gradient));
	if ( (abs(gradient[0]) < epsilon) and (abs(gradient[1]) < epsilon) and (abs(gradient[2]) < epsilon)):
	# if ( ((x- epsilon) < 0) and ((y - epsilon) < 0) and ((z - epsilon) < 0)):
		return [0,0,0];

	# return gradient;
	return helper.normalizeVector(gradient);
	# return helper.normalizeVector([x,y,z]);

def markImage(img, _i, _j):
	i = int(_i);
	j = int(_j);
	# print(i,j);
	# print(len(img), len(img[0]))
	img[i][j] = 255;
	img[i][j][1] = 0;
	img[i][j][2] = 0;
	return img;

rootPath = "/Users/joelpp/Documents/Maitrise/Code/experiments/";
expNumber = int(sys.argv[1]);
path, dirs, files = os.walk(rootPath+repr(expNumber)+"/probes/").next();
output = open(rootPath+repr(expNumber)+"/gradientList.txt", 'wb');


# for index in xrange(1):
# for index in xrange(602,603):
for index in xrange(len(files)):
	img = misc.imread(rootPath+repr(expNumber)+"/probes/Probe_"+repr(index)+".png");
	# img = misc.imread(rootPath+repr(expNumber)+"/renders/Render_"+repr(index)+".png");

	height, width = helper.getImgSize(img);

	gradient = getGradient(img);

	print("Obtained gradient for probe "+repr(index)+". grad = "+repr(gradient));
	output.write(repr(gradient[0])+" "+repr(gradient[1])+" "+repr(gradient[2])+"\n");

	(phi, theta) = helper.cartesianToSpherical_yup(gradient[0], gradient[1], gradient[2]);
	(u,v) = helper.equiPTtoUV_yup(phi, theta);
	(i,j) = helper.uvToIJ_yup(u,v, width, height);
	# helper.show("phi", phi);
	# helper.show("theta", theta);
	# helper.show("u", u);
	# helper.show("v", v);
	# helper.show("i", i);
	# helper.show("j", j);
	img = markImage(img,i,j);
	misc.imsave(rootPath+"/"+repr(expNumber)+"/misc/Probe_"+repr(index)+".png", img);