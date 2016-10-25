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
    WeightMatrixType generateWeightsMatrix(int NumberOfSamples);
	Eigen::VectorXd generatebVector(int NumberOfSamples, String optimizationFolderPath);
	void outputWeightsMatrixToFile(int NumberOfSamples, String optimizationFolderPath);
    void outputBVectorToFile(int NumberOfSamples, String optimizationFolderPath);

	float generateTriplet(int row, int col, const ProbeInterpolationRecord* iRec);
	void generateTriplets(int NumberOfSamples, String outputPath, std::vector<Eigen::Triplet<float>>* eigenTriplets, bool optimizeForCoeffs);
	void generateRGBValuesFromProbes(int NumberOfSamples, String savePath, Eigen::VectorXd* eigenVector);
    void generateRGBValuesFromProbes(int NumberOfSamples);
    void generateRGBValuesFromSamples(int NumberOfSamples, String savePath, Eigen::VectorXd* eigenVector = 0);
    void generateInterpolatedCoefficientsFromProbes(int NumberOfSamples, String savePath, Eigen::VectorXd* eigenVector = 0);
	void save();
	void load(int maxSamples);
	std::vector<float> tryOptimizationPass(int NumberOfSamples, bool optimizeForMitsubaSamples, String optimizationFolderPath);
	bool probeOptimizationPass(WeightMatrixType& A, Eigen::VectorXd& b, Eigen::VectorXd* result);
    void createbVector(Eigen::VectorXd* bVector, const Eigen::VectorXd* rgbColumn, String& optimizationFolderPath);

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