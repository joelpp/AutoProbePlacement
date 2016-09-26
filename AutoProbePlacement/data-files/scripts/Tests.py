import helper;

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def printTestResult(result):
	if (result == False):
		print bcolors.FAIL + "Test failed :(" + bcolors.ENDC;
	if (result == True):
		print bcolors.OKGREEN + "Test passed :)" + bcolors.ENDC;
# IJ to UV conversion should convert screen-space (i,j) coordinates to the range [-1,1]


def testIJtoUV():
	print "Testing helper.ijToUV_yup"
	width = 256;
	height = 512;
	
	(i,j) = (0,0);
	(u,v) = helper.ijToUV_yup(i, j, height, width);
	
	print "(i,j):" + repr((i,j)),
	print "-> (u,v):" + repr((u,v))
	printTestResult((u,v) == (-1,-1));
	print;
	(i,j) = (512,256);
	(u,v) = helper.ijToUV_yup(i, j, height, width);
	print "(i,j):" + repr((i,j)),
	print "-> (u,v):" + repr((u,v))
	printTestResult((u,v) == (1,1));
	
def testUVToPT():
	print "Testing helper.equiUVtoPT_yup"

	
	(u,v) = (0,0);
	(phi,theta) = helper.equiUVtoPT_yup(u,v);
	print "(u,v):" + repr((u,v)) + "-> (phi,theta):" + repr((phi,theta))
	printTestResult((phi,theta) == (0, helper.pi/2.0));
	print;

	(u,v) = (-1,0);
	(phi,theta) = helper.equiUVtoPT_yup(u,v);
	print "(u,v):" + repr((u,v)) + "-> (phi,theta):" + repr((phi,theta))
	printTestResult((phi,theta) == (-helper.pi, helper.pi/2.0));
	print;
	
	(u,v) = (0.5,0);
	(phi,theta) = helper.equiUVtoPT_yup(u,v);
	print "(u,v):" + repr((u,v)) + "-> (phi,theta):" + repr((phi,theta))
	printTestResult((phi,theta) == (helper.pi/2, helper.pi/2.0));
	print;
	
	(u,v) = (1,0);
	(phi,theta) = helper.equiUVtoPT_yup(u,v);
	print "(u,v):" + repr((u,v)) + "-> (phi,theta):" + repr((phi,theta))
	printTestResult((phi,theta) == (helper.pi, helper.pi/2.0));
	print;
	
	(u,v) = (0,-1);
	(phi,theta) = helper.equiUVtoPT_yup(u,v);
	print "(u,v):" + repr((u,v)) + "-> (phi,theta):" + repr((phi,theta))
	printTestResult((phi,theta) == (0,0));
	print;
	
def testSphericalToCartesian():
	print "Testing helper.sphericalToCartesian_yup"
	(r,phi,theta) = (1,0,0);
	(x,y,z) = helper.sphericalToCartesian_yup(r,phi,theta);
	print "(r,phi,theta):" + repr((r,phi,theta)) + "-> (x,y,z):" + repr((x,y,z))
	printTestResult((x,y,z) == (0,1,0));
	print;
	
	(r,phi,theta) = (1,0,helper.pi / 2.0);
	(x,y,z) = helper.sphericalToCartesian_yup(r,phi,theta);
	print "(r,phi,theta):" + repr((r,phi,theta)) + "-> (x,y,z):" + repr((x,y,z))
	printTestResult((abs(x-1) < 1e-5) and (abs(y-0) < 1e-5) and (abs(z-0) < 1e-5));
	print;
	
	(r,phi,theta) = (1,helper.pi/2,helper.pi / 2.0);
	(x,y,z) = helper.sphericalToCartesian_yup(r,phi,theta);
	print "(r,phi,theta):" + repr((r,phi,theta)) + "-> (x,y,z):" + repr((x,y,z))
	printTestResult((abs(x-0) < 1e-5) and (abs(y-0) < 1e-5) and (abs(z-1) < 1e-5));
	print;
	
def testCartesianToSpherical():
	print "Testing helper.cartesianToSpherical_yup"
	(x,y,z) = (0,1,0);
	(phi, theta) = helper.cartesianToSpherical_yup(x,y,z);
	print "(x,y,z):" + repr((x,y,z)) + "-> (phi, theta):" + repr((phi, theta))
	printTestResult((abs(phi-0) < 1e-5) and (abs(theta-0) < 1e-5));
	print;

	(x,y,z) = (1,0,0);
	(phi, theta) = helper.cartesianToSpherical_yup(x,y,z);
	print "(x,y,z):" + repr((x,y,z)) + "-> (phi, theta):" + repr((phi, theta))
	printTestResult((abs(phi-0) < 1e-5) and (abs(theta-helper.pi/2) < 1e-5));
	print;
	
	(x,y,z) = (0,0,1);
	(phi, theta) = helper.cartesianToSpherical_yup(x,y,z);
	print "(x,y,z):" + repr((x,y,z)) + "-> (phi, theta):" + repr((phi, theta))
	printTestResult((abs(phi-helper.pi/2) < 1e-5) and (abs(theta-helper.pi/2) < 1e-5));
	print;
	
def testSHxyz_yup():
	print "Testing helper.SHxyz_yup"
	(x,y,z) = (1,0,0);
	(L,M) = (0,0);
	SH = helper.SHxyz_yup((L,M),(x,y,z));
	print "((L,M), x,y,z):" + repr(((L,M),(x,y,z))) + "-> (SH):" + repr(SH)
	printTestResult(abs(SH-0.282095) < 1e-5);
	print;
	
	(L,M) = (1,-1);
	SH = helper.SHxyz_yup((L,M),(x,y,z));
	print "((L,M), x,y,z):" + repr(((L,M),(x,y,z))) + "-> (SH):" + repr(SH)
	printTestResult(abs(SH-0) < 1e-5);
	print;
	
	(L,M) = (1,0);
	SH = helper.SHxyz_yup((L,M),(x,y,z));
	print "((L,M), x,y,z):" + repr(((L,M),(x,y,z))) + "-> (SH):" + repr(SH)
	printTestResult(abs(SH-0) < 1e-5);
	print;
	
	(L,M) = (1,1);
	SH = helper.SHxyz_yup((L,M),(x,y,z));
	print "((L,M), x,y,z):" + repr(((L,M),(x,y,z))) + "-> (SH):" + repr(SH)
	printTestResult(abs(SH-0) < 1e-5);
	print;
###########################	
	
print
# testIJtoUV();
# testUVToPT();
# testSphericalToCartesian();
testCartesianToSpherical();
# testSHxyz_yup();