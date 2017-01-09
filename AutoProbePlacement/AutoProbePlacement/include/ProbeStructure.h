#include <G3D/G3DAll.h>
#include <unordered_map>
#include <map>

#ifndef PROBESTRUCTURE_H
#define PROBESTRUCTURE_H

#include "Tetrahedron.h"
#include "Probe.h"
#include <fstream>
#include "Actor.h"

struct SProbe
{
    float position[4];
    float coefficients[27];
    float gradients[81];
};

struct SProbeStructure
{
    SProbe probes[343];
    int dimensions[4]; // int3 padded into int4
    float firstProbePosition[4];  // vec3 padded into vec4
    float step;
};


struct ProbeInterpolationRecord
{
	G3D::Array<int> probeIndices;
	G3D::Array<float> weights;
};

enum EProbeStructureType
{
	Trilinear = 0,
	Closest = 1,
	Tetrahedral = 2,
	WeightedNearestNeighbour = 3,
	NUM_TYPES = 4
};

namespace std {
    template <>
        class hash<Vector3>{
        public :
            size_t operator()(const Vector3 &vec ) const
            {
                return hash<int>()( (int) vec.x ) ^ hash<int>()( (int)vec.y ) ^ hash<int>()( (int)vec.z );
            }
    };
};

class ProbeStructure {
protected:
	
	Array<LineSegment> lineArray;
	Array<Tetrahedron*> tetrahedraArray;
	Array<Array<int> > tetFaces;
	Array<Vector3> tetFaceNormals;
	Array<Vector3> tetFaceNormalPositions;
	String probeStructurePath;

public:
	
	static Array<G3D::String> typeMap;
	EProbeStructureType m_type;
	inline void createTypeMap();

	std::unordered_map<Vector3, Probe*> probeMap;

	ProbeStructure();
	// TODO : remove path argument and build it in the constructor instead sheesh
	ProbeStructure(String sceneName, String probeStructureName);
	ProbeStructure(String sceneName, String probeStructureName, int numProbes, EProbeStructureType type);
	void operator=(const ProbeStructure &p );
	void reset();
	void makeLineArray();
	void makeTetrahedraArray();
	void loadHullFaces();
	void computeTetFaceNormals();
	void computeTetVertexNormals();
	void makeProbeList();
	void loadTetrahedraNeighbors();
	bool writeProbePositionsToFile();
	bool writeProbePositionsToTetgenFile();
	Vector3 getFaceNormalPosition(int i);
	Vector3 getFaceNormal(int i);
	int tetFaceNormalCount();
	LineSegment getLineSegment(int i);
	int lineArraySize();
	Tetrahedron* getTetrahedron(int i);
	Probe* getProbe(int i);
	int getProbeIndex(const G3D::Vector3& coord);
	int probeCount();
	void changeProbeCoeff(int probeIndex, int coeffIndex, Vector3 amount);
	void changeProbePosition(int probeIndex, Vector3 dp);
	void updateProbes(bool updateAll);
	void applyOffsetToProbes(std::vector<float>& displacements);
	bool isOutsideSceneBounds(G3D::Vector3 pos, float tolerance);
	void displaceProbesWithGradient(std::vector<float>& displacements, float maxStepLength);


	String name() { return m_name; }

	String type() { return String(typeMap[m_type]); }

	EProbeStructureType eType() { return m_type; }

	float gamma() { return m_gamma; }

    void setGamma(float gamma);

    int width() { return m_width; }

    void setWidth(int w) { m_width = w; }

    int height() { return m_height; }

    void setHeight(int h) { m_height = h; }

    int numSamples() { return m_NumSamples; }

    void setNumSamples(int n) { m_NumSamples = n; }

	float step() { return m_step;  }

	void setStep(float step);

	//void coefficientInterpolation(G3D::Vector3 position, std::vector<G3D::Vector3> interpolatedCoefficients, std::vector<int>& probeIndices, std::vector<float> weights);
	ProbeInterpolationRecord getInterpolationProbeIndicesAndWeights(const G3D::Vector3& position);
	G3D::Array<G3D::Vector3> ProbeStructure::getInterpolatingProbesCoords(const G3D::Vector3& pos, int step);
	G3D::Array<int> getInterpolatingProbeIndices(const G3D::Vector3& pos);


	std::fstream probeListFileHandle(bool reading);
	std::fstream infoFileHandle(bool reading);
    std::fstream probeCoefficientsFileHandle(int i, bool reading);
    std::fstream probeCoefficientsGradientsFileHandle(int i, bool reading);

	TProbeCoefficients interpolatedCoefficients(const G3D::Vector3& position, const G3D::Vector3& normal, int NumberOfCoeffs);

	Array<G3D::Vector3> ProbeStructure::reconstructSHPerBand(const G3D::Vector3& position, const G3D::Vector3& normal, int NumberOfCoeffs);

    G3D::Vector3 reconstructSH(const G3D::Vector3& position, const G3D::Vector3& normal, int NumberOfCoeffs);

	void addProbe(const G3D::Vector3& position);

	void generateProbes(std::string type, bool allProbes, bool generateGradients, bool showOutput);

    void extractSHCoeffs(bool generateGradients, bool bUploadToGPU);

	void loadProbeStructureInfo();

	void loadSceneInfo();

    void savePositions(bool useManipulator);

    void uploadToGPU();

	void setIntegrator(String integrator);

	void setType(String type);

	void saveInfoFile();

    bool hasProbes();

	void saveCoefficients();

	void updateAll(bool showOutput);

	void deleteAllProbes();

	void removeProbe(int i);

	void createDirectoryTree();

	std::vector<int> m_dimensions;
	Array<Probe*> probeList;
	String m_name;
	String m_sceneName;
	String m_integrator;

    float m_step;
    float m_firstProbePosition[3];
    float m_gamma;
	float m_gradientDisplacement;

    int m_NumCoefficients;
    int m_NumColors;
    int m_width;
    int m_height;
    int m_NumSamples;


};

#endif

