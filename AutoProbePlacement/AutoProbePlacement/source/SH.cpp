#include "SH.h"

#define PI 3.141592654f

#define AXIS_X 0
#define AXIS_Y 1
#define AXIS_Z 2

int SH::factorial(int n)
{
	if (n == 0)
	{
		return 1;
	}
	int fact = n;
	for (int i = n - 1; i >= 1; --i)
	{
		fact *= i;
	}

	return fact;
}

float SH::S(int m, float x, float y)
{
	float** M = new float*[m + 1];
	for (int i = 0; i < m + 1; i++)
		M[i] = new float[2];

	//float M[m + 1][2];
	M[0][0] = 0;
	M[0][1] = 1;

	for (int i = 1; i < m + 1; ++i)
	{
		M[i][0] = x * M[i - 1][0] + y * M[i - 1][1];
		M[i][1] = x * M[i - 1][1] - y * M[i - 1][0];
	}
	float toReturn = M[m][0];

	for (int i = 0; i < m + 1; i++)
	{
		delete[] M[i];
	}
	delete[] M;

	return toReturn;
}

float SH::C(int m, float x, float y)
{
	float** M = new float*[m + 1];
	for (int i = 0; i < m + 1; i++)
		M[i] = new float[2];

	//float M[m + 1][2];
	M[0][0] = 0;
	M[0][1] = 1;

	for (int i = 1; i < m + 1; ++i)
	{
		M[i][0] = x * M[i - 1][0] + y * M[i - 1][1];
		M[i][1] = x * M[i - 1][1] - y * M[i - 1][0];
	}
	float toReturn = M[m][1];

	for (int i = 0; i < m + 1; i++)
	{
		delete[] M[i];
	}
	delete[] M;

	return toReturn;
}

float SH::K(int l, int m)
{
	return sqrt((2 * l + 1) * factorial(l - abs(m)) / (4 * PI * factorial(l + abs(m))));
}

float SH::legendreP(int lMax, int mMax, float z)
{
	float** P = new float*[lMax + 1];
	for (int i = 0; i < lMax + 1; i++)
		P[i] = new float[mMax + 1];
	P[0][0] = 1;

	for (int m = 0; m < mMax + 1; ++m)
	{
		for (int l = m; l < lMax + 1; ++l)
		{
			if (m == l)
			{
				if (m == 0) //#if we have P00, continue
				{
					continue;
				}
				//#otherwise, calculate Pmm
				P[m][m] = (1 - 2 * m) * P[m - 1][m - 1];
			}
			else if (l == m + 1)
			{
				P[l][m] = (2 * m + 1) * z * P[m][m];
			}
			else
			{
				P[l][m] = ((2 * l - 1) * z * P[l - 1][m] - (l + m - 1) * P[l - 2][m]) / (l - m);
			}
		}
	}
	float toReturn = P[lMax][mMax];

	for (int i = 0; i < lMax + 1; i++)
	{
		delete[] P[i];
	}
	delete[] P;

	return toReturn;
}

float SH::SHxyz_yup(int l, int m, const G3D::Vector3& normal)
{
	const float& x = normal[0];
	const float& y = normal[1];
	const float& z = normal[2];

	if (m == 0)
	{
		return K(l, m) * legendreP(l, m, y);
	}
	else if (m < 0)
	{
		return sqrtf(2.f) * K(l, m) * S(-m, x, z) * legendreP(l, -m, y);
	}
	else
	{
		return sqrtf(2.f) * K(l, m) * C(m, x, z) * legendreP(l, m, y);
	}
}


std::pair<int, int> SH::kToLM(int k)
{
	std::pair<int, int> lm;

	lm.first = (int) floor(sqrtf((float) k));
	int offset = k - (int) powf((float) lm.first, 2.f);
	lm.second = -lm.first + offset;

	return lm;
}

SHGradient SH::gradient(int NumCoefficients)
{
	SHGradient grad;
	for (int i = 0; i < 3; ++i)
	{
		grad.append(G3D::Array<float>());

		for (int j = 0 ; j < NumCoefficients; ++j)
		{
			grad[i].append(0);
		}
	}

	return grad;

}

SHGradient SH::gradients(Vector3 dir)
{
	int NumCoefficients = 9;

	SHGradient grad = gradient(NumCoefficients);

	float x = dir[0];
	float y = dir[1];
	float z = dir[2];

	float x2pz2 = pow(x,2) + pow(z,2);
    float EPSILON = 1e-5;

	for (int i =0 ; i < NumCoefficients; ++i)
	{
		std::pair<int, int> lm = kToLM(i);
		int l = lm.first;
		int m = lm.second;
		float sqrt1Plus2l = sqrtf(1 + 2 * l);
		float sqrtPi = sqrt(PI);
		float sqrt2Pi = sqrt(2 * PI);
		float cosMaXZ = cos( m * atan2(x, z) );
		float sinMaXZ = sin( m * atan2(x, z) );
		float factLminusAbsM = (float) factorial(l - abs(m));
		float factLplusAbsM = (float) factorial(l + abs(m));
		
		if (m == 0)
		{

			float legendrePLy = legendreP(l, 0, y);
			float legendreP1Ly = legendreP(l + 1, 0, y);

			if (x2pz2 > EPSILON)
			{
				grad[AXIS_X][i] = -( (1 + l) * sqrt1Plus2l * x * y * (y * legendrePLy - legendreP1Ly) ) / (2.f * sqrtPi * x2pz2);
				grad[AXIS_Z][i] = -( (1 + l) * sqrt1Plus2l * y * z * (y * legendrePLy - legendreP1Ly) ) / (2.f * sqrtPi * x2pz2);
			}
			else
			{
				grad[AXIS_X][i] = 0;
				grad[AXIS_Z][i] = 0;
			}
			
			grad[AXIS_Y][i] = ( (1 + l) * sqrt1Plus2l * (y * legendrePLy - legendreP1Ly) ) / (2.f * sqrtPi);
		}
		else if (m > 0)
		{
			float legendreP1Lmy = legendreP(l + 1, m, y);
			float legendrePlmy = legendreP(l, m, y);
			
			if (x2pz2 > EPSILON)
			{
				grad[AXIS_X][i] = ( (sqrt1Plus2l * factLminusAbsM) / factLplusAbsM ) * ( (1 + l - m) * x * y * cosMaXZ * legendreP1Lmy + legendrePlmy * ( -((1 + l) * x * pow(y, 2) * cosMaXZ) + m * z * sinMaXZ) ) / (sqrt2Pi * x2pz2);
				grad[AXIS_Z][i] = -( (sqrt1Plus2l * factLminusAbsM) / factLplusAbsM) * ( (-1 - l + m) * y * z * cosMaXZ * legendreP1Lmy +  legendrePlmy * ( (1 + l) * z *pow(y, 2) * cosMaXZ + m * x  * sinMaXZ) ) / (sqrt2Pi * x2pz2);
			}
			else
			{
				grad[AXIS_X][i] = 0;
				grad[AXIS_Z][i] = 0;
			}	
			grad[AXIS_Y][i] = ((cosMaXZ * sqrt1Plus2l * factLminusAbsM) / factLplusAbsM) * ((1 + l) * y * legendrePlmy + (-1 - l + m) * legendreP1Lmy) / (sqrt2Pi);
		}
		else if (m < 0)
		{
			float legendreP1Lminusmy = legendreP(l + 1, -m, y);
			float legendrePlminusmy = legendreP(l, -m, y);

			if (x2pz2 > EPSILON)
			{
				grad[AXIS_X][i] = ( (sqrt1Plus2l * factLminusAbsM) / factLplusAbsM ) * ( -((1 + l + m) * x * y * legendreP1Lminusmy * sinMaXZ ) + legendrePlminusmy * ( m * z  * cosMaXZ + (1 + l) * x *pow(y, 2) * sinMaXZ) ) / (sqrt2Pi * x2pz2);
				grad[AXIS_Z][i] = ( (sqrt1Plus2l * factLminusAbsM) / factLplusAbsM ) * ( -( (1 + l + m) * y * z * legendreP1Lminusmy * sinMaXZ ) + legendrePlminusmy * (-(m * x  * cosMaXZ) + (1 + l) * pow(y, 2) * z * sinMaXZ)) / (sqrt2Pi * x2pz2);
			}
			else
			{
				grad[AXIS_X][i] = 0;
				grad[AXIS_Z][i] = 0;
			}	
			grad[AXIS_Y][i] = ( (sqrt1Plus2l * factLminusAbsM) / factLplusAbsM ) * ( -( (1 + l) * y * legendrePlminusmy ) + (1 + l + m) * legendreP1Lminusmy) * sinMaXZ / sqrt2Pi;
		}
	}
	return grad;
}
