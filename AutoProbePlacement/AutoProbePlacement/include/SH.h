
#ifndef SH_H
#define SH_H

#include <G3D/G3DAll.h>
#include <math.h>
typedef Array<Array<float>> SHGradient;


class SH{
protected:
	static SHGradient gradient(int NumCoefficients);
		
public:
	static int factorial(int n);
	static float S(int m, float x, float y);
	static float C(int m, float x, float y);
	static float K(int l, int m);
	static float legendreP(int lMax, int mMax, float z);
	static float SHxyz_yup(int l, int m, const G3D::Vector3& normal);
	static std::pair<int, int> kToLM(int k);
	static SHGradient gradients(Vector3 dir);
};


#endif