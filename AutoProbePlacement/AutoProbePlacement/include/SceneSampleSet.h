#include <vector>
#include "Eigen/Sparse"
#include "Eigen/SparseQR"
#include <G3D/G3DAll.h>
#include <GLG3D/GLG3D.h>
typedef Eigen::SparseMatrix<double> WeightMatrixType;

struct ProbeInterpolationRecord;
class SceneSample;
class ProbeStructure;

enum ESSFile
{
    Samples = 0,
    Values = 1
};

class SceneSampleSet
{
public:
	static bool create(String sceneName, String sampleSetName);

	SceneSampleSet();
	SceneSampleSet(std::string sceneName, std::string sampleSetName, float scale, int numSamplesToLoad);
	
	void addSample(SceneSample sample);
	void generateTriplets(int NumberOfSamples, std::vector<Eigen::Triplet<float>>* eigenTriplets);
	void generateRGBValuesFromProbes(int NumberOfSamples, bool ref, Eigen::VectorXd* eigenVector);
    void generateRGBValuesFromProbes(int NumberOfSamples);
    void save();
	void load(int maxSamples);
	std::vector<float> tryOptimizationPass(int NumberOfSamples, bool ref);
	bool probeOptimizationPass(WeightMatrixType& A, Eigen::VectorXd& b, Eigen::VectorXd* result);
	void createbVector(Eigen::VectorXd* bVector, const Eigen::VectorXd* rgbColumn);

	void clearValues();
	void clearPositions();

    std::fstream SceneSampleSet::openFile(ESSFile type, bool reading);

	/*
		Member variables
	*/
	G3D::Array<G3D::Point3> m_points;
	G3D::Array<G3D::Color3> m_colors;
	std::vector<SceneSample> m_samples;
	std::string m_sceneName, m_sampleSetName;
	float m_scale;
	ProbeStructure* probeStructure;
	bool outputToLog;

};