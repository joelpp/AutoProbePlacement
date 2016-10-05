#include <G3D/G3DAll.h>
#include <unordered_map>
#include <map>

#ifndef PROBESTRUCTURE_H
#define PROBESTRUCTURE_H

#include "Tetrahedron.h"
#include "Probe.h"
#include <fstream>
#include "Actor.h"

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
	
	static std::vector<std::string> typeMap;
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
	void displaceProbesWithGradient(std::vector<float>& displacements);


	String name() { return m_name; }

	String type() { return String(typeMap[m_type]); }

	float gamma() { return m_gamma; }

	//void coefficientInterpolation(G3D::Vector3 position, std::vector<G3D::Vector3> interpolatedCoefficients, std::vector<int>& probeIndices, std::vector<float> weights);
	ProbeInterpolationRecord getInterpolationProbeIndicesAndWeights(G3D::Vector3 position);
	G3D::Array<G3D::Vector3> ProbeStructure::getInterpolatingProbesCoords(const G3D::Vector3& pos, int step);
	G3D::Array<int> getInterpolatingProbeIndices(const G3D::Vector3& pos);

	void loadProbeInfo();

	float m_step;
	std::vector<int> m_dimensions;
	float m_firstProbePosition[3];
	Array<Probe*> probeList;
	String m_name;
	String m_sceneName;
	float m_gamma;

};

#endif

