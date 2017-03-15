import math
# import numpy
pi = 3.141592654

def computeSHGradients(dir):
	dYdx = [0.0] * 9;
	dYdy = [0.0] * 9;
	dYdz = [0.0] * 9;
	x = dir[0];
	y = dir[1];
	z = dir[2];
	x2pz2 = (x ** 2 + z ** 2);
	epsilon = 1e-5;

	for i in xrange(9):
		(l,m) = kToLM_new(i);
		sqrt1Plus2l = math.sqrt(1 + 2 * l);
		sqrtPi = math.sqrt(pi);
		sqrt2Pi = math.sqrt(2 * pi);
		cosMaXZ = math.cos( m * math.atan2(x, z) )
		sinMaXZ = math.sin( m * math.atan2(x, z) )
		factLminusAbsM = math.factorial(l - math.fabs(m))
		factLplusAbsM = math.factorial(l + math.fabs(m))
		
		if (m == 0):
			legendrePLy = legendreP(l, 0, y);
			legendreP1Ly = legendreP(l + 1, 0, y);

			if ((x**2 + z** 2) > epsilon):
				dYdx[i] = -( (1 + l) * sqrt1Plus2l * x * y * (y * legendrePLy - legendreP1Ly) ) / (2. * sqrtPi * x2pz2);
				dYdz[i] = -( (1 + l) * sqrt1Plus2l * y * z * (y * legendrePLy - legendreP1Ly) ) / (2. * sqrtPi * x2pz2);
			else:
				dYdx[i] = 0;
				dYdz[i] = 0;		
			
			dYdy[i] = ( (1 + l) * sqrt1Plus2l * (y * legendrePLy - legendreP1Ly) ) / (2. * sqrtPi);

		elif (m > 0):
			legendreP1Lmy = legendreP(l + 1, m, y);
			legendrePlmy = legendreP(l, m, y);
			
			if ((x**2 + z** 2) > epsilon):
				dYdx[i] = ( (sqrt1Plus2l * factLminusAbsM) / factLplusAbsM ) * ( (1 + l - m) * x * y * cosMaXZ * legendreP1Lmy + legendrePlmy * ( -((1 + l) * x * y ** 2 * cosMaXZ) + m * z * sinMaXZ) ) / (sqrt2Pi * x2pz2);
				dYdz[i] = -( (sqrt1Plus2l * factLminusAbsM) / factLplusAbsM) * ( (-1 - l + m) * y * z * cosMaXZ * legendreP1Lmy +  legendrePlmy * ( (1 + l) * z * y ** 2 * cosMaXZ + m * x  * sinMaXZ) ) / (sqrt2Pi * x2pz2);
			else:
				dYdx[i] = 0;
				dYdz[i] = 0;
				
			dYdy[i] = ((cosMaXZ * sqrt1Plus2l * factLminusAbsM) / factLplusAbsM) * ((1 + l) * y * legendrePlmy + (-1 - l + m) * legendreP1Lmy) / (sqrt2Pi);

		elif (m < 0):
			legendreP1Lminusmy = legendreP(l + 1, -m, y);
			legendrePlminusmy = legendreP(l, -m, y);

			if ((x**2 + z** 2) > epsilon):
				dYdx[i] = ( (sqrt1Plus2l * factLminusAbsM) / factLplusAbsM ) * ( -((1 + l + m) * x * y * legendreP1Lminusmy * sinMaXZ ) + legendrePlminusmy * ( m * z  * cosMaXZ + (1 + l) * x * y ** 2 * sinMaXZ) ) / (sqrt2Pi * x2pz2);
				dYdz[i] = ( (sqrt1Plus2l * factLminusAbsM) / factLplusAbsM ) * ( -( (1 + l + m) * y * z * legendreP1Lminusmy * sinMaXZ ) + legendrePlminusmy * (-(m * x  * cosMaXZ) + (1 + l) * y ** 2 * z * sinMaXZ)) / (sqrt2Pi * x2pz2)
			else:
				dYdx[i] = 0;
				dYdz[i] = 0;
				
			dYdy[i] = ( (sqrt1Plus2l * factLminusAbsM) / factLplusAbsM ) * ( -( (1 + l) * y * legendrePlminusmy ) + (1 + l + m) * legendreP1Lminusmy) * sinMaXZ / sqrt2Pi;
	# dYdx = computedYdx(dir);
	# dYdy = computedYdy(dir);
	# dYdz = computedYdz(dir);
	return (dYdx, dYdy, dYdz);
	
def computedYdx( dir ):

	dYdx = [0.0] * 9;
	x = dir[0];
	y = dir[1];
	z = dir[2];
	
	t1 = x*x;
	t2 = y*y;
	t3 = z*z;
	t4 = t1+t2+t3;
	t5 = math.sqrt(t4);
	t7 = 1.0/t5/t4;
	t19 = 1.0/t4;
	t23 = t4*t4;
	t24 = 1.0/t23;
	
	dYdx[0] = 0.0;
	dYdx[1] = -0.488603*y*t7*x;
	dYdx[2] = -0.488603*z*t7*x;
	dYdx[3] = 0.488603/t5-0.488603*t1*t7;
	dYdx[4] = 0.1092548E1*y*t19-0.2185096E1*t1*y*t24;
	dYdx[5] = -0.2185096E1*x*y*t24*z;
	dYdx[6] = -0.1892352E1*t3*t24*x;
	dYdx[7] = 0.1092548E1*z*t19-0.2185096E1*z*t1*t24;
	dYdx[8] = 0.1092548E1*x*t19-0.1092548E1*(t1-t2)*t24*x;
	return dYdx;

def computedYdy( dir ):
	dYdy = [0.0] * 9;
	x = dir[0];
	y = dir[1];
	z = dir[2];

	t1 = x*x;
	t2 = y*y;
	t3 = z*z;
	t4 = t1+t2+t3;
	t5 = math.sqrt(t4);
	t9 = 1/t5/t4;
	t13 = y*t9;
	t18 = 1/t4;
	t22 = t4*t4;
	t23 = 1/t22;
	
	dYdy[0] = 0.0;
	dYdy[1] = 0.488603/t5-0.488603*t2*t9;
	dYdy[2] = -0.488603*t13*z;
	dYdy[3] = -0.488603*t13*x;
	dYdy[4] = 0.1092548E1*x*t18-0.2185096E1*x*t2*t23;
	dYdy[5] = 0.1092548E1*z*t18-0.2185096E1*t2*z*t23;
	dYdy[6] = -0.1892352E1*y*t3*t23;
	dYdy[7] = -0.2185096E1*x*y*t23*z;
	dYdy[8] = -0.1092548E1*y*t18-0.1092548E1*(t1-t2)*t23*y;
	
	return dYdy;

def computedYdz( dir ):
	dYdz = [0.0] * 9;
	x = dir[0];
	y = dir[1];
	z = dir[2];

	t1 = x*x;
	t2 = y*y;
	t3 = z*z;
	t4 = t1+t2+t3;
	t5 = math.sqrt(t4);
	t7 = 1/t5/t4;
	t20 = t4*t4;
	t21 = 1/t20;
	t25 = 1/t4;
	
	dYdz[0] = 0.0;
	dYdz[1] = -0.488603*y*t7*z;
	dYdz[2] = 0.488603/t5-0.488603*t3*t7;
	dYdz[3] = -0.488603*z*t7*x;
	dYdz[4] = -0.2185096E1*x*y*t21*z;
	dYdz[5] = 0.1092548E1*y*t25-0.2185096E1*y*t3*t21;
	dYdz[6] = 0.1892352E1*z*t25-0.1892352E1*t3*z*t21;
	dYdz[7] = 0.1092548E1*x*t25-0.2185096E1*t3*t21*x;
	dYdz[8] = -0.1092548E1*(t1-t2)*t21*z;
	
	return dYdz;


#Get values for the first 9 SH bands (from PP Sloan)
def SH(L, M, x, y, z):
	invpi = 1 / math.sqrt(math.pi);
	if (L == 0):
		return invpi/(2);
	elif (L == 1):
		if (M == -1):
			return -invpi * math.sqrt(3) * y / 2;

		elif (M == 0):
			return invpi * math.sqrt(3) * z / 2;
		elif (M == 1):
			return -invpi * math.sqrt(3) * x / 2;
	elif (L == 2):
		if (M == -2):
			return invpi * math.sqrt(15) * x * y / 2;

		elif (M == -1):
			return -invpi * math.sqrt(15) * z * y / 2;
		elif (M == 0):
			return invpi * math.sqrt(5) * (3*z*z -1) / 4;
		elif (M == 1):
			return -invpi * math.sqrt(15) * z * x / 2;
		elif (M == 2):
			return invpi * math.sqrt(15) * (x*x - y*y) / 4;
	elif (L == 3):
		if (M == -3):
			return invpi * math.sqrt(2*35) * (3*x*x - y*y) * y / 8;
		elif (M == -2):
			return invpi * math.sqrt(105) * x * y * z / 2;
		elif (M == -1):
			return -invpi * math.sqrt(2*21) * (-1 + 5*z*z) * y / 8;
		elif (M == 0):
			return invpi * math.sqrt(7) * z * (5*z*z - 3) / 4;
		elif (M == 1):
			return -invpi * math.sqrt(2*21) * (5*z*z - 1) * x / 8;
		elif (M == 2):
			return invpi * math.sqrt(105) * (x*x - y*y) * z / 4;
		elif (M == 3):
			return -invpi * math.sqrt(2*35) * (-3*y*y + x*x) * x / 8;
	elif (L == 4):
		if (M == -4):
			return 3*invpi * math.sqrt(35) * (x*x - y*y) * y * x / 4;		
		elif (M == -3):
			return -3 * invpi * math.sqrt(2*35) * (3*x*x - y*y) * y * z / 8;
		elif (M == -2):
			return 3 * invpi * math.sqrt(5) * x * y * (-1+7*z*z) / 4;
		elif (M == -1):
			return -3*invpi * math.sqrt(2*5) * (-3 + 7*z*z) * y * z / 8;
		elif (M == 0):
			return 3 * invpi * (35*z*z*z*z - 30*z*z + 3) / 16;
		elif (M == 1):
			return -3 * invpi * math.sqrt(2*5) * (-3 + 7*z*z) * x * z / 8;
		elif (M == 2):
			return 3 * invpi * math.sqrt(5) * (x*x - y*y) * (-1+7*z*z) / 8;
		elif (M == 3):
			return -3 * invpi * math.sqrt(2*35) * (-3*y*y + x*x) * x * z / 8;
		elif (M == 4):
			return 3 * invpi * math.sqrt(35) * (x*x*x*x - 6*y*y*x*x + y*y*y*y) / 16;

def SHlookup(i,j,coeff):
	return shArray[i][j][coeff];

def SHBetterLookup(i,j,coeff):
	# print("Trying to grab coeff="+repr(coeff));
	if (shArray[coeff] == None):
		# shArray[coeff] = cPickle.load(open('/Users/joelpp/Documents/Maitrise/Code/IrradianceMapping/SHvalues/values_'+repr(coeff)+'.p', 'rb'));
		# shArray[coeff] = cPickle.load(open('/Users/joelpp/Documents/Maitrise/Code/IrradianceMapping/SHvalues_yup/512x256/values_'+repr(coeff)+'.p', 'rb'));
		shArray[coeff] = cPickle.load(open('../SH/128x64/values_'+repr(coeff)+'.p', 'rb'));
	
	return shArray[coeff][i][j];

def SHBetterLookup_small(i,j,coeff):
	# print("Trying to grab coeff="+repr(coeff));
	if (shArray[coeff] == None):
		# shArray[coeff] = cPickle.load(open('/Users/joelpp/Documents/Maitrise/Code/IrradianceMapping/SHvalues/values_'+repr(coeff)+'.p', 'rb'));
		# shArray[coeff] = cPickle.load(open('/Users/joelpp/Documents/Maitrise/Code/IrradianceMapping/SHvalues_yup/256x128/values_'+repr(coeff)+'.p', 'rb'));

		shArray[coeff] = cPickle.load(open('../SH/128x64/values_'+repr(coeff)+'.p', 'rb'));
	# print "i "+repr(i)
	# print "j "+repr(j)
	# print "coeff "+repr(coeff)
	return shArray[coeff][i][j];

def SHxyz(L,M, x,y,z):
	if (M == 0):
		return K(L,M) * legendreP(L,M,z);
	elif (M < 0):
		return math.sqrt(2) * K(L,M) * S(-M, x,y) * legendreP(L,-M, z);
	else:
		return math.sqrt(2) * K(L,M) * C(M, x,y) * legendreP(L,M, z);

# def SHxyz_yup((L,M), x,y,z):
	# if (M == 0):
		# return K(L,M) * legendreP(L,M,y);
	# elif (M < 0):
		# return math.sqrt(2) * K(L,M) * S(-M, x,z) * legendreP(L,-M, y);
	# else:
		# return math.sqrt(2) * K(L,M) * C(M, x,z) * legendreP(L,M, y);

def SHxyz_yup(L,M, normal):
	x = normal[0];
	y = normal[1];
	z = normal[2];
	if (M == 0):
		return K(L,M) * legendreP(L,M,y);
	elif (M < 0):
		return math.sqrt(2) * K(L,M) * S(-M, x,z) * legendreP(L,-M, y);
	else:
		return math.sqrt(2) * K(L,M) * C(M, x,z) * legendreP(L,M, y);

def S(m,x,y):
	returnArray = numpy.zeros((m+1,2));
	returnArray[0][0] = 0;
	returnArray[0][1] = 1;

	for i in range(1,m+1):
		returnArray[i][0] = x * returnArray[i-1][0] + y * returnArray[i-1][1];
		returnArray[i][1] = x * returnArray[i-1][1] - y * returnArray[i-1][0];

	return returnArray[m][0];

def C(m,x,y):
	returnArray = numpy.zeros((m+1,2));
	returnArray[0][0] = 0;
	returnArray[0][1] = 1;

	for i in range(1,m+1):
		returnArray[i][0] = x * returnArray[i-1][0] + y * returnArray[i-1][1];
		returnArray[i][1] = x * returnArray[i-1][1] - y * returnArray[i-1][0];

	return returnArray[m][1];



def legendreP(l_,m_,z):
	P = numpy.zeros((l_+1, m_+1));
	# print(P);
	P[0][0] = 1;

	for m in range(m_+1):
		# show("m",m);
		for l in range(m,l_+1):
			# show("l",l);
			if (m == l):
				if (m == 0): #if we have P00, continue
					continue;
				#otherwise, calculate Pmm
				P[m][m] = (1 - 2*m) * P[m-1][m-1];
			elif (l == m+1):
				P[l][m] = (2*m + 1) * z * P[m][m];
			else:
				P[l][m] = ((2*l - 1) * z * P[l-1][m] - (l+m-1)*P[l-2][m]) / (l-m);
			# show("P", P[l][m]);


	return P[l_][m_];



def SHRam(L, M, x,y,z):
	invpi = 1 / math.sqrt(math.pi);
	if (L == 0):
		return 0.282095;
	elif (L == 1):
		if (M == -1):
			return 0.488603 * y;

		elif (M == 0):
			return 0.488603 * z;
		elif (M == 1):
			return 0.488603 * x;
	elif (L == 2):
		if (M == -2):
			return 1.092548 * x * y;

		elif (M == -1):
			return 1.092548 * y * z;
		elif (M == 0):
			return 0.315392 * (3*z*z - 1);
		elif (M == 1):
			return 1.092548 * x * z;
		elif (M == 2):
			return 0.546274 * (x*x - y*y);


def lmToK(l,m):
	if (l == 0):
		return 0;
	elif (l == 1):
		return m+2;
	elif (l == 2):
		return m + 6;
def kToLM_new(k):
	l = math.floor(math.sqrt(k));
	offset = k - l*l;
	return (int(l), int(-l + offset));

#Get a tuple representing the SH bands from an index.
def kToLM(k):
	if (k == 0):
		return (0,0);
	elif (k<4):
		return (1, k-2);
	elif (k < 9):
		return (2, k-6);
	elif (k < 16):
		return (3, k - 12);
	elif (k < 25):
		return (4, k - 20);
	elif (k < 36):
		return (5, k - 30);
	elif (k < 49):
		return (6, k - 42);
	elif (k < 64):
		return (7, k - 56);
	elif (k < 81):
		return (8, k - 72);
	elif (k < 100):
		return (9, k - 90);
	elif (k < 121):
		return (10, k - 110);

#Get a tuple representing the SH bands from an index.
def kToLMLegendre(k):
	if (k == 0):
		return (0,0);
	elif (k<3):
		return (1, k-1);
	elif (k < 6):
		return (2, k-3);
	elif (k < 10):
		return (3, k-6);
	else:
		return (4, k - 20);
def K(l, m):
	# if ((l - math.fabs(m)) < 0):
		# print(l,m);
	return math.sqrt( (2*l + 1) * math.factorial(l - math.fabs(m)) / (4*math.pi * math.factorial(l+ math.fabs(m))) );
