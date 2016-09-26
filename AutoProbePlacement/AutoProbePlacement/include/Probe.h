#ifndef PROBE_H
#define PROBE_H

#include <G3D/G3DAll.h>

class ProbeManipulator;

typedef Array<Array<float>> SHGradient;
typedef Array<Array<Vector3>> CoeffGradients;



enum EResource
{
	Probes,
	Positions,
	Normals,
	Coefficients,
	CoefficientGradients
};
class Probe {
protected:


public:
	static shared_ptr<ArticulatedModel> model;

	Probe();
	Probe(int i, String probeStructurePath);
	shared_ptr<Texture> getTexture(int texID);
	void computeCoefficientsFromTexture(bool alsoSet);
	String buildPath(EResource res);

	shared_ptr<Texture>		texture;
	shared_ptr<Texture>		irr_texture;
	Point3					position;
	Array<Vector3>			coeffs;
	Array<Array<Vector3> >	coeffGradients;
	Vector3					gradient;
	Vector3					normal;
	Sphere					m_sphere;
	String probeStructurePath;
	int index;

	shared_ptr<ProbeManipulator> manipulator;
	shared_ptr<ProbeManipulator> getManipulator();
	Point3 getPosition();
	CFrame frame;


	bool bNeedsUpdate;
};
class ProbeManipulator : public ThirdPersonManipulator
{
public:
	Probe* p;
	static shared_ptr<ProbeManipulator> ProbeManipulator::create() {
		return shared_ptr<ProbeManipulator>(new ProbeManipulator());
	}

	virtual void onDrag(const Vector2 &delta)
	{
		p->bNeedsUpdate = true;
		ThirdPersonManipulator::onDrag(delta);
	}

};
#endif