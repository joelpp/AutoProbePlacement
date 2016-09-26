from scipy import misc
import numpy
import cPickle
import math
import time;
import os;
import sys
from helper import *

#------------------------------
#|		USEFUL FUNCTIONS	  |	  	
#------------------------------

def interpolatedCoeffs():
	return [[0.22835, 0.175914, 0.0521269], [-0.166509, 0.0396026, -0.0546848], [0.00143171, 0.00135783, 0.000510839], [0.150878, 0.0976164, -0.0590761], [-0.0741768, 0.215812, 0.0792114], [-0.00117166, 2.08309e-05, -0.000672084], [-0.109057, -0.114924, -0.0530471], [0.000138019, 0.000225896, -0.000776082], [0.10825, 0.0941369, 0.0167043], [-0.185625, 0.0407479, -0.0701099], [0.000635143, 0.00206264, 0.00117922], [0.0452764, 0.0213999, 0.0444015], [-0.00092091, -0.00130549, -0.00102507], [0.0521527, 0.0286284, 0.0530418], [0.000501433, 0.00054247, 0.000245857], [0.0326711, 0.0473007, 0.0400749], [-0.0339273, 0.0782164, 0.0181387], [-0.00139628, -0.000605961, -0.00110602], [-0.0866144, -0.0651651, -0.061446], [0.00074802, 0.00100933, 0.00103525], [0.0294984, 0.039697, 0.045034], [0.00179317, 0.00158871, 0.0012901], [0.00698199, 0.000471491, -0.0120407], [0.000475711, 0.000427855, 0.000606637], [-0.101219, -0.0798149, -0.0649765], [0.00817279, 0.0543438, 0.0242529], [0.000323149, 0.000450322, 0.000370197], [0.0127613, 0.0664099, 0.0460535], [-0.00208395, -0.00126988, -0.00163712], [-0.0302183, -0.0415847, -0.0354316], [0.0011959, 0.000960789, 0.00132436], [-0.0524038, -0.0606827, -0.0446377], [-0.000116663, 0.000127417, -0.00036197], [-0.0207234, -0.0103462, -0.0248949], [-0.00121615, -0.0012125, -0.00106536], [0.0326314, 0.0446376, 0.0469873], [-0.0560336, -0.0113127, -0.0394774], [0.000402749, 0.000379321, 0.0003714], [-0.0366481, 0.0102621, -0.0152859], [0.00106124, 0.00191738, 0.00134648], [0.041469, 0.0357625, 0.0465079], [-0.00128434, -0.00134142, -0.00116014], [-0.0401214, -0.0290482, -0.0357306], [-0.00127895, -0.00179808, -0.00150282], [0.0168824, 0.00440556, 0.0110969], [-0.000724738, -0.000446141, -0.000652565], [0.0292138, 0.0393272, 0.03798], [0.000805309, 0.000888835, 0.000817694], [0.00117487, 0.00158159, -0.00840869], [-0.00847576, 0.0524573, 0.0196431], [-0.000625083, -0.000468918, -0.000613217], [-0.01971, 0.00607013, -0.0104428], [-0.000577457, -0.000134498, -0.00052075], [-0.0481214, -0.0379476, -0.0337797], [0.00131615, 0.00153325, 0.00165613], [0.0297481, 0.0276591, 0.0256393], [-0.00153668, -0.00134422, -0.00136372], [0.0281788, 0.0342758, 0.0338223], [0.00067484, 0.000529287, 0.000468232], [0.019082, 0.0184039, 0.0148688], [0.000919364, 0.00099214, 0.00106186], [-0.032168, -0.0311978, -0.0262881], [-8.75001e-05, -4.98774e-05, -0.000152057], [-0.0107538, -0.0138719, -0.0211326], [-0.0201246, 0.0369553, 0.00845112], [0.000185295, 0.000419096, 0.000242633], [0.00655263, 0.0244239, 0.0148986], [-0.000144505, -5.49757e-05, -0.00017823], [-0.000165983, 0.0177213, 0.0125298], [-0.00158353, -0.00116293, -0.00130043], [-0.0301943, -0.0322723, -0.0313437], [0.0010431, 0.00103049, 0.00102704], [0.0246757, 0.0263646, 0.0252879], [0.00147716, 0.00133187, 0.00143232], [-0.00829685, -0.0135709, -0.0100215], [0.000442411, 0.000646477, 0.000449762], [-0.0240971, -0.017277, -0.0216645], [-0.000782217, -0.000800943, -0.000746218], [0.00797323, 0.00472102, 0.00551028], [-0.000310312, -0.000381335, -0.000344172], [0.0215033, 0.022795, 0.0231682], [-0.0441287, -0.00573905, -0.0253038], [0.000243704, 0.000349881, 0.000292087], [-0.017693, 0.0057449, -0.00475422], [0.000260641, 0.000337066, 0.000263216], [-0.00636457, 0.00844431, 0.00042869], [0.000511574, 0.000890108, 0.000593863], [0.0162089, 0.0192697, 0.0221599], [-0.00151966, -0.00137336, -0.00133028], [-0.0139363, -0.0157845, -0.0156365], [0.00102442, 0.00114718, 0.0011476], [-0.0265918, -0.0221727, -0.0230781], [-0.000376166, -0.000612471, -0.00053108], [-0.0020577, -0.00732624, -0.00547892], [-0.000785285, -0.000626129, -0.000697337], [0.0104786, 0.0136092, 0.0137534], [0.000223354, 0.000255611, 0.000169884], [0.0100826, 0.00972666, 0.00598865], [0.000282166, 0.000191513, 0.000299413], [-0.0110414, -0.0048256, -0.00368752]];

#Get the coefficients for Grace Cathedral's irradiance maps from Ramamoorthi
def loadGraceCoeff():
	# return [[0.078908,  0.043710,  0.054161],[0.039499,  0.034989,  0.060488],[0.033974, -0.018236, -0.026940],[0.029213, -0.005562,  0.000944],[0.011141, -0.005090, -0.012231],[0.026240, -0.022401, -0.047479],[0.015570, -0.009471, -0.014733],[0.056014,  0.021444,  0.013915],[0.021205, -0.005432, -0.030374]];
	return [[0.78908,  0.43710,  0.54161],[0.39499,  0.34989,  0.60488],[-0.33974, -0.18236, -0.26940],[-0.29213, -0.05562,  0.00944],[-0.11141, -0.05090, -0.12231],[-0.26240, -0.22401, -0.47479],[-0.15570, -0.09471, -0.14733],[0.56014,  0.21444,  0.13915],[0.21205, -0.05432, -0.30374]];

def loadEucalyptusCoeff():
	return [[0.38,0.43,0.45],[0.29,0.36,0.41],[0.04,0.03,0.01],[-.10,-.10,-.09],[-0.06,-0.06,-0.04],[0.01,-0.01,-0.05],[-0.09,-0.13,-0.15],[-0.06,-0.05,-0.04],[0.02,0,-0.05]];

def diffuseCoeffs((l,m)):
	if (l == 0):
		return 3.14159;
	elif (l == 1):
		return 2.094395;
	elif (l == 2):
		return 0.785398;

def phongCoeffs((l,m), r):
	pi = 3.14159;
	if (l == 0):
		return pi;
	elif (l == 1):
		return pi * (1+r)/(2+r);
	elif (l == 2):
		return pi * r / (3+r);
	elif (l == 3):
		return pi * (r-1)*(r+1) / (8 + 6*r + r*r);
	elif (l == 4):
		return (r-2)*r / ((r+3)*(r+5));
	elif (l == 5):
		return (r-3) * (r-1) * (r+1) / ((2+r)*(4+r)*(6+r));
	elif (l == 6):
		return (r-4) * (r-2) * (r) / ((3+r)*(5+r)*(7+r));
	elif (l == 7):
		return (r-5) * (r-3) * (r-1) * (r+1) / ((2+r)*(4+r)*(6+r)*(8+r));
	elif (l == 8):
		return (r-6) * (r-4) * (r-2) * r / ((3+r)*(5+r)*(7+r)*(9+r));
	elif (l == 9):
		return (r-7) * (r-5) * (r-3) * (r-1) * (r+1) / ((2+r)*(4+r)*(6+r)*(8+r)*(10+r));

	# elif (l == 10):
	# 	return 0.03125 * (r-8) * (r-6) * (r-4) / ((2+r)*(4+r)*(6+r));
#------------------------------
#|		PROGRAM START 		  |	  	
#------------------------------
start = time.time();
# path = sys.argv[1];
# print(path);
# print(sphericalToCartesian(1,0,0));
# img = misc.imread("/Users/joelpp/Documents/Maitrise/Images/"+path)
# img = misc.imread("/Users/joelpp/Documents/Maitrise/Images/Mitsuba/SphericalCamera_BetterScaleObj_0_3.png")
# height = len(img);
# width = len(img[0]);
# height = 128;
# width = 256;
height = 64;
width = 128;
maxBand = 2;
numCoeffs = 0;
nullList = [];
for i in range(maxBand+1):
	numCoeffs = numCoeffs + (2*i + 1);

for i in xrange(numCoeffs):
	nullList.append([0,0,0]);
phongExponent = 1;
expNumber = int(sys.argv[1]);
rootPath = "/Users/joelpp/Documents/Maitrise/Code/experiments/";

path, dirs, files = os.walk(rootPath+repr(expNumber)+"/coefficients/").next();

# probeNumber = 156;
# Load a custom SH coefficients file
# numx = int(sys.argv[1]);
# numy = int(sys.argv[2]);
# numz = int(sys.argv[3]);
offset = (0,0,0);
# for s in range(numx):
# 	for yindex in range(numy):
# 		for t in range(numz):
		# print("s="+repr(s)+" , v="+repr(yindex)+", t="+repr(t));
for probeNumber in xrange(len(files)):
	savePath = "/Users/joelpp/Documents/Maitrise/code/experiments/"+repr(expNumber)+"/renders/Render_"+repr(probeNumber)+".png";
	if (os.path.isfile(savePath)):
		print("Probe #"+repr(probeNumber)+" already rendered! continuing...");
		continue;

	L = cPickle.load(open('/Users/joelpp/Documents/Maitrise/Code/experiments/'+repr(expNumber)+'/coefficients/Coefficients_'+repr(probeNumber)+'.p', 'rb'));
	# L = cPickle.load(open('/Users/joelpp/Documents/Maitrise/Images/grace/coefficients.p', 'rb'));
	# L = interpolatedCoeffs();

	# LMC = cPickle.load(open('MonteCarloCoefficientsSloan.p', 'rb'));

	# LRam = cPickle.load(open('QuadratureCoefficientsRamamoorthi_'+repr(s)+'_'+repr(t)+'.p', 'rb'));
	# LMCRam = cPickle.load(open('MonteCarloCoefficientsRamamoorthi.p', 'rb'));

	# for i in range(9):
	# 	for j in range(3):
	# 		L[i][j] *= 1;
	# 		# LMC[i][j] *= 1;
	# 		# LMCRam[i][j] *= 1;

	# L = loadGraceCoeff();
	# LRam = loadGraceCoeff();
	# L = loadEucalyptusCoeff();
	print("Probe #"+repr(probeNumber));
	print(type(L));
	print(L);
	# print(nullList);


	# Weights for the SH bands, as described in Ramamoorthi's article
	ramaCoeffs = numpy.zeros((5));
	ramaCoeffs[0] = 0.429043;
	ramaCoeffs[1] = 0.511664;
	ramaCoeffs[2] = 0.743125;
	ramaCoeffs[3] = 0.886227;
	ramaCoeffs[4] = 0.247708;

	numChannels = 3;
	#Initialize image arrays.
	RGB = numpy.zeros((height, width, numChannels));
	# RGBmc = numpy.zeros((height, width, 3));
	# RGBRam = numpy.zeros((height, width, 3));
	# RGBMcRam = numpy.zeros((height, width, 3));
	# thetas = numpy.zeros((height, width,3));
	# phis = numpy.zeros((height, width,3));
	# xs = numpy.zeros((height, width));
	# ys = numpy.zeros((height, width));
	# zs = numpy.zeros((height, width));
	# tester = numpy.zeros((height, width));
	# SHpix = numpy.zeros((9, height, width,3));

	if (L == nullList):
		print("Coeffs are empty. Continuing...");
		misc.imsave(savePath, RGB);
		continue;

	#Initialize some values for the counter.
	total = height*width;
	perten = math.floor(total / 10);
	progress = 0;
	epsilon = 0.0001;

	# Loop over all them texels.
	for i in range(height):
		for j in range(width):

			#Increment progress counter if appropriate
			if ((i*width+j) % perten == 0):
				print(repr(progress)+"% complete")
				progress += 10;

			#Convert from i,j to u,v coordinates
			# (u,v) = ijToUV(i,j, height, width);
			(u,v) = ijToUV_yup(i,j, height, width);

			#Convert from u,v to theta,phi coordinates
			# (phi, theta) = equiUVtoPT(u,v);
			(phi, theta) = equiUVtoPT_yup(u,v);

			#Convert from theta,phi to x,y,z coordinates
			# (x,y,z) = sphericalToCartesian(1,phi, theta);
			(x,y,z) = sphericalToCartesian_yup(1,phi, theta);

			#For validation purposes: draw the normalized value of each independant cartesian coordinate to its own array
			# xs[i,j] = (x+1)/2;
			# ys[i,j] = (y+1)/2;
			# # zs[i,j] = (z+1)/2;

			# print("SH(1,1): "+repr(SH((2,1), x, y, z)));
			# print("SHpt(1,1): "+repr(SHpt((2,1), phi, theta)));

			# Custom test for more validation.
			# tester[i,j] = 1.0;
			# if ((phi == math.pi/2.0) or (theta == 0)):
			# 	# show("theta", theta);
			# 	tester[i,j] = 0.0;

			# Validation: write to their own array the normalized values of theta and phi.
			# for k in range(3):
			# 	thetas[i,j,k] = theta * 255 / (2*math.pi)
			# 	phis[i,j,k] = phi * 255 / (1*math.pi)

			# Validation: write an RGB value corresponding to each SH band's value at that point; blue is more negative and red is more positive
			# for k in range(9):
			# 	value = SH(kToLM(k), x, y, z);
			# 	SHpix[k,i,j,0] = (value + 1) / 2;
			# 	SHpix[k,i,j,2] = 1.0 - SHpix[k,i,j,0];

			# Color the point according to the SH coefficients as described in Ramamoorthi's article.

			if (numChannels == 4):
				RGB[i][j][3] = 1.0;
			for k in range(3):
				# This is the rendering equation using Sloan's coefficients
				# RGB[i][j][k] = 0.4288 * L[8][k] * (x*x - y*y) + 0.247 * 3 * L[6][k] * z*z + 2.7834 * L[0][k] - 0.247 * L[6][k] + 2.6943 * ( L[4][k]*x*y - L[7][k]*x*z - L[5][k]*y*z) + 3.2149 * (-L[3][k] * x - L[1][k]*y + L[2][k]*z);
				
				for coeff in range(numCoeffs):
					# RGB[i][j][k] += L[coeff][k] * diffuseCoeffs(kToLM(coeff)) * SH(kToLM(coeff),x,y,z);
					# if (type(SHlookup(i,j,coeff) is not float):
					# 	show("Wrong type!", type(SHlookup(i,j,k)));
					# 	show("i", i);
					# 	show("j", j);
					# 	show("k", k)
					# # THE GOOD PART
					# RGB[i][j][k] += L[coeff][k] * phongCoeffs(kToLM(coeff), 1) * SHxyz(kToLM(coeff), x,y,z);
					RGB[i][j][k] += L[coeff][k] * phongCoeffs(kToLM(coeff), phongExponent) * SHBetterLookup(i,j,coeff);

				# RGBmc[i][j][k] = 0.4288 * LMC[8][k] * (x*x - y*y) + 0.247 * 3 * LMC[6][k] * z*z + 2.7834 * LMC[0][k] - 0.247 * LMC[6][k] + 2.6943 * ( LMC[4][k]*x*y - LMC[7][k]*x*z - LMC[5][k]*y*z) + 3.2149 * (-LMC[3][k] * x - LMC[1][k]*y + LMC[2][k]*z);
				
				# This is the rendering equation using Ramamoorthi's coefficients
				# RGBRam[i][j][k] = ramaCoeffs[0] * LRam[8][k] * (x*x - y*y) + ramaCoeffs[2] * LRam[6][k] * z*z + ramaCoeffs[3] * LRam[0][k] - ramaCoeffs[4] * LRam[6][k] + 2.0 * ramaCoeffs[0] * ( LRam[4][k]*x*y + LRam[7][k]*x*z + LRam[5][k]*y*z) + 2.0 * ramaCoeffs[1] * (LRam[3][k] * x + LRam[1][k]*y + LRam[2][k]*z);
				# RGBMcRam[i][j][k] = ramaCoeffs[0] * LMCRam[8][k] * (x*x - y*y) + ramaCoeffs[2] * LMCRam[6][k] * z*z + ramaCoeffs[3] * LMCRam[0][k] - ramaCoeffs[4] * LMCRam[6][k] + 2.0 * ramaCoeffs[0] * ( LMCRam[4][k]*x*y + LMCRam[7][k]*x*z + LMCRam[5][k]*y*z) + 2.0 * ramaCoeffs[1] * (LMCRam[3][k] * x + LMCRam[1][k]*y + LMCRam[2][k]*z);

	# Write arrays into images.
	# misc.imsave("/Users/joelpp/Documents/Maitrise/Images/IrradianceMaps/MonteCarloRenderSloan.png", RGBmc);

	misc.imsave(savePath, RGB);
	# misc.imsave("/Users/joelpp/Documents/Maitrise/Images/grace/Render.png", RGB);
		
	# misc.imsave("/Users/joelpp/Documents/Maitrise/Images/IrradianceMaps/MonteCarloRenderRamamoorthi.png", RGBMcRam);
	# misc.imsave("/Users/joelpp/Documents/Maitrise/Images/IrradianceMaps/QuadratureRenderRamamoorthi_"+repr(s)+"_"+repr(t)+".png", RGB);
# misc.imsave("images/phis.png", phis);
# misc.imsave("images/thetas.png", thetas);
# misc.imsave("images/xs.png", xs);
# misc.imsave("images/ys.png", ys);
# misc.imsave("images/zs.png", zs);
# misc.imsave("images/tester.png", tester);

# for k in range(9):
# 	misc.imsave("images/SH_"+repr(kToLM(k))+".png", SHpix[k]);

# knock knock. Time to go!
end = time.time();
print("Elapsed: "+repr((end-start)));
