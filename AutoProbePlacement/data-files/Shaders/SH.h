#define PI 3.141592654


int[2] kToLM(int k)
{
	int toReturn[2];

	int l = int(floor(sqrt(k)));
	int offset = k - l*l;
	toReturn[0] = l;
	toReturn[1] = -l + offset;

    return toReturn;
}

int factorial(int n)
{
	if (n == 0)
	{
		return 1;
	}
	int fact = n;
	int i = n;
	while (i >= 1)
	{
		fact *= i;
		i--;
	}

	return fact;
}
float S(int m, float x, float y)
{
	float M[10][2];
	M[0][0] = 0;
	M[0][1] = 1;

	for (int i = 1; i < m + 1; ++i)
	{
		M[i][0] = x * M[i - 1][0] + y * M[i - 1][1];
		M[i][1] = x * M[i - 1][1] - y * M[i - 1][0];
	}

	return M[m][0];
}

float C(int m, float x, float y)
{
	float M[10][2];
	M[0][1] = 1;
	M[0][0] = 0;

	for (int i = 1; i < m + 1; ++i)
	{
		M[i][0] = x * M[i - 1][0] + y * M[i - 1][1];
		M[i][1] = x * M[i - 1][1] - y * M[i - 1][0];
	}

	return M[m][1];
}

float K(int l, int m)
{
	return sqrt((2 * l + 1) * factorial(l - abs(m)) / (4 * PI * factorial(l + abs(m))));
}

float legendreP(int lMax, int mMax, float z)
{
	float P[10][10];
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

	return P[lMax][mMax];
}

float SHxyz_yup(int l, int m, vec3 normal)
{
	float x = normal[0];
	float y = normal[1];
	float z = normal[2];

	if (m == 0)
	{
		return K(l, m) * legendreP(l, m, y);
	}
	else if (m < 0)
	{
		return sqrt(2) * K(l, m) * S(-m, x, z) * legendreP(l, -m, y);
	}
	else
	{
		return sqrt(2) * K(l, m) * C(m, x, z) * legendreP(l, m, y);
	}
}

float SHFunctionValue(uint shBand, vec3 normal)
{
	int[2] lm = kToLM(int(shBand));
	
	return SHxyz_yup(lm[0], lm[1], normal);
}