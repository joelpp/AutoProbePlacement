
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
	static double S(int m, double x, double y);
	static double C(int m, double x, double y);
	static double K(int l, int m);
	static double legendreP(int lMax, int mMax, double z);
	static double SHxyz_yup(int l, int m, const G3D::Vector3& normal);
	static std::pair<int, int> kToLM(int k);
	static SHGradient gradients(Vector3 dir);
};


#endif