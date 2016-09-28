#include <vector>
#include "Eigen/Sparse"
#include "Eigen/SparseQR"
typedef Eigen::SparseMatrix<double> WeightMatrixType;
namespace G3D
{
	class Vector3;
	class ArticulatedModel;
	template <class T, size_t MIN_ELEMENTS = 10>
	class Array;
}

struct ProbeInterpolationRecord;
class SceneSample;
class ProbeStructure;



class SceneSampleSet
{
public:
	std::vector<SceneSample> m_samples;
	std::string m_sceneName, m_sampleSetName;
	float m_scale;
	ProbeStructure* probeStructure;
	bool outputToLog;

	SceneSampleSet();	
	SceneSampleSet(std::string sceneName, std::string sampleSetName, float scale);
	
	void addSample(SceneSample sample);
	void generateTriplets(int NumberOfSamples, std::vector<Eigen::Triplet<float>>* eigenTriplets);
	void generateRGBValuesFromProbes(int NumberOfSamples, bool ref, Eigen::VectorXd* eigenVector);
	void save();
	void load();
	std::vector<float> tryOptimizationPass(int NumberOfSamples, bool ref);
	bool probeOptimizationPass(WeightMatrixType& A, Eigen::VectorXd& b, Eigen::VectorXd* result);
	void createbVector(Eigen::VectorXd* bVector, const Eigen::VectorXd* rgbColumn);

	void clear();
};