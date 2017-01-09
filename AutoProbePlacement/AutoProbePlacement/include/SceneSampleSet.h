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
    Values = 1,
	Coeffs
};

class SceneSampleSet
{
public:
	static bool create(String sceneName, String sampleSetName);

	SceneSampleSet();
	SceneSampleSet(std::string sceneName, std::string sampleSetName, float scale, int numSamplesToLoad);
	
	void addSample(SceneSample sample);
    WeightMatrixType generateWeightsMatrix(int NumberOfSamples, int NumberOfCoeffs);
	Eigen::VectorXd generatebVector(int NumberOfSamples, int NumberOfCoeffs, String optimizationFolderPath);
	void outputWeightsMatrixToFile(int NumberOfSamples, int NumberOfCoeffs, String optimizationFolderPath);
    void outputBVectorToFile(int NumberOfSamples, int NumberOfCoeffs, String optimizationFolderPath);

	float generateTriplet(int row, int col, const ProbeInterpolationRecord* iRec);
	void generateTriplets(int NumberOfSamples, int NumberOfCoeffs, String outputPath, std::vector<Eigen::Triplet<float>>* eigenTriplets, bool optimizeForCoeffs);
	void generateRGBValuesFromProbes(int NumberOfSamples, int NumberOfCoeffs, String savePath, Eigen::VectorXd* eigenVector);
    void generateRGBValuesFromProbes(int NumberOfSamples, int NumberOfCoeffs);
    void generateRGBValuesFromSamples(int NumberOfSamples, String savePath, Eigen::VectorXd* eigenVector = 0);
	void generateInterpolatedCoefficientsFromProbes(int NumberOfSamples, int NumberOfCoeffs);
	void generateInterpolatedCoefficientsFromProbes(int NumberOfSamples, int NumberOfCoeffs, String savePath, Eigen::VectorXd* eigenVector = 0);
	void save();
	bool load(int maxSamples);
	std::vector<float> tryOptimizationPass(int NumberOfSamples, int NumberOfCoeffs, bool optimizeForMitsubaSamples, String optimizationFolderPath);
	bool probeOptimizationPass(WeightMatrixType& A, Eigen::VectorXd& b, Eigen::VectorXd* result);
    void createbVector(Eigen::VectorXd* bVector, const Eigen::VectorXd* rgbColumn, String& optimizationFolderPath);

	void clearValues();
	void clearPositions();

	float distanceToProbe(const G3D::Vector3& position, int p);
	float w(const G3D::Vector3& position, const G3D::Vector3& normal, int m, int n, int axis, int color);
	float A(const G3D::Vector3& position, const G3D::Vector3& normal, int m, int n, int axis, int color);
	float dWeightMdProbeN(const G3D::Vector3& position, const G3D::Vector3& normal, int m, int n, int axis, int color);
	float dRdX(const G3D::Vector3& position, const G3D::Vector3& normal, int NumberOfCoeffs, int m, int axis, int color);
	float R(const G3D::Vector3& position, const G3D::Vector3& normal, int NumberOfCoeffs, int m, int axis, int color);
	float B(const G3D::Vector3& position, const G3D::Vector3& normal, int NumberOfCoeffs, int m, int n, int axis, int color);
	float dInverseSquaredSumdProbeN(const G3D::Vector3& position, const G3D::Vector3& normal, int m, int n, int axis, int color);

	float D(const G3D::Vector3& position, const G3D::Vector3& normal, int m, int n, int axis, int color);
	float InverseSumOf1OverSquaredProbeDistances(const G3D::Vector3& position);
	float SumOfInverseSquaredProbeDistances(const G3D::Vector3& position);
	float dInverseDistanceSquaredMdProbeN(const G3D::Vector3& position, const G3D::Vector3& normal, int m, int n, int axis, int color);
    std::fstream openFile(ESSFile type, bool reading);
	void removeDarkSamples();
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
	bool coeffReference; //as opposed to RGB reference
	bool oneRowPerSHBand; //as opposed to RGB reference

};