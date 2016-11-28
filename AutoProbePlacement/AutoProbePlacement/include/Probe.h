#ifndef PROBE_H
#define PROBE_H

#include <G3D/G3DAll.h>

class ProbeManipulator;

typedef Array<Array<float>> SHGradient;
typedef Array<Array<Vector3>> CoeffGradients;
typedef Array<Vector3> TProbeCoefficients;

enum EResource
{
	Probes,
	Probes_X_NEG,
	Probes_X_POS,
	Probes_Y_NEG,
	Probes_Y_POS,
	Probes_Z_NEG,
	Probes_Z_POS,
	Positions,
	Normals,
	Coefficients,
	CoefficientGradients,

};
class Probe {
protected:


public:
	static shared_ptr<ArticulatedModel> model;

	Probe();
	Probe(int i, String probeStructurePath);
	shared_ptr<Texture> getTexture(int texID);
	void reconstructSH(const G3D::Vector3& normal);
	void computeCoefficientsFromTexture(bool alsoSet, bool computeGradients, float gradientDisplacement);
	String buildPath(EResource res);
	void saveCoefficients();
	void computeCoefficients(std::shared_ptr<G3D::Image> probeTexture, TProbeCoefficients& coeffs);
	CoeffGradients computeProbeCoeffGradients(float gradientDisplacement);

	shared_ptr<Texture>		texture;
	shared_ptr<Texture>		irr_texture;
	Point3					position;
	TProbeCoefficients		coeffs;
	Array<Array<Vector3> >	coeffGradients;
	Vector3					gradient;
	Vector3					normal;
	Sphere					m_sphere;
	String probeStructurePath;
	int index;

	shared_ptr<ProbeManipulator> manipulator;
	shared_ptr<ProbeManipulator> getManipulator();
	Point3 getPosition();
    void setPosition(G3D::Vector3& pos);
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