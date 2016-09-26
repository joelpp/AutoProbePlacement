import numpy as np
from scipy.misc import *


a = np.load("./result.npy");

img = np.zeros((len(a), len(a[0]), 3));
print type(img);
r = 1;
g = 1; 
b = 1;

for j in xrange(len(img)):
	for i in xrange(len(img[0])):
		img[i][j][0] = a[i][j] * r;
		img[i][j][1] = a[i][j] * g;
		img[i][j][2] = a[i][j] * b;
imsave('./outfile.png', img)


# print a;