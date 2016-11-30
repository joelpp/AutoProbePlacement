#include <iostream>
#include <fstream>
//#include <utility>
#include <string>
#include <sstream>

#include "SceneSampleSet.h"
#include "SceneSample.h"
#include "ProbeStructure.h"
#include "SH.h"
#include "Helpers.h"


#define COLOR_RED	0
#define COLOR_GREEN 1
#define COLOR_BLUE	2

#define AXIS_X 0
#define AXIS_Y 1
#define AXIS_Z 2


#define WEIGHTS_TRILERP 0
#define WEIGHTS_CLOSEST 1
#define WEIGHTS_TETRAHEDRAL 0

#define USE_SPARSE 1

#define GENERATE_POSITION_FITTING_MATRIX 1
#define GENERATE_COEFFICIENTS_FITTING_MATRIX 0

//#define USE_CUSTOM_FUNCTION
#define ONE_ROW_PER_SH_BAND

bool invertNormal = true;
void probeOptimizationPass(std::vector<int> matrixRows, std::vector<int> matrixColumns, std::vector<float> matrixElements, std::vector<float> rgb);


SceneSampleSet::SceneSampleSet(std::string sceneName, std::string sampleSetName, float scale, int numSamplesToLoad)
{
	m_sceneName = sceneName;
	m_sampleSetName = sampleSetName;
	m_scale = scale;
	if (!load(numSamplesToLoad))
	{
		throw std::exception("Sample set load failed");
	}
}
bool SceneSampleSet::load(int maxSamples)
{
	m_samples.clear();
	std::fstream samplesFile, irradianceFile;
	std::string samplesPath = "../data-files/Scenes/" + m_sceneName + "/SampleSets/" + m_sampleSetName + "/SamplePositions.txt";

	bool loadIrradianceResults2 = true;
	std::string irradiancePath;

	if (loadIrradianceResults2)
	{
		irradiancePath = "../data-files/Scenes/" + m_sceneName + "/SampleSets/" + m_sampleSetName + "/IrradianceResults2.txt";
	}
	else
	{
		irradiancePath = "../data-files/Scenes/" + m_sceneName + "/SampleSets/" + m_sampleSetName + "/IrradianceResults.txt";

	}

	samplesFile.open(samplesPath.c_str(), std::fstream::in);
	irradianceFile.open(irradiancePath.c_str(), std::fstream::in);

	if (!samplesFile.is_open())
	{
		debugPrintf("Couldn't open position samples file!\n");
		return false;
	}

	if (!irradianceFile.is_open())
	{
		debugPrintf("Couldn't open irradiance results file!\n");
		return false;
	}

	int NumberOfLinesPerSample = 6;
	int currentLine = 0;
	std::string sampleLine, irradianceLine;
	SceneSample ss;
	int numLoaded = 0;
	while (std::getline(samplesFile, sampleLine))
	{
		Array<String> splitLine = stringSplit(String(sampleLine.c_str()), ' ');

		if (currentLine % NumberOfLinesPerSample == 0)
		{
			ss = SceneSample();

			if (!std::getline(irradianceFile, irradianceLine))
			{
                ss.irradiance = Color3::zero();
            }
            else
            {
                Array<String> splitLine = stringSplit(String(irradianceLine.c_str()), ' ');

                float divider = 1.f;
                ss.irradiance = Color3(std::stof(splitLine[0].c_str()) / divider,
                                       std::stof(splitLine[1].c_str()) / divider,
                                       std::stof(splitLine[2].c_str()) / divider);
            }
			
		}

		if (currentLine == 2)
		{ // Samples Positions
			ss.position = Vector3(std::stof(splitLine[0].c_str()) * m_scale,
								  std::stof(splitLine[1].c_str()) * m_scale,
								  std::stof(splitLine[2].c_str()) * m_scale);
		}

		if (currentLine == 4)
		{ // Normals
			ss.normal = Vector3(std::stof(splitLine[0].c_str()),
				std::stof(splitLine[1].c_str()),
				std::stof(splitLine[2].c_str()));
		}

		currentLine++;

		if (currentLine == NumberOfLinesPerSample)
		{
			m_samples.push_back(ss);
			numLoaded++;
			if (numLoaded == maxSamples)
			{
				break;
			}
			currentLine = 0;
		}

	}

	samplesFile.close();
	irradianceFile.close();

	m_points.clear();
	m_colors.clear();
	for (SceneSample& ss : m_samples)
	{
		if (ss.irradiance != Color3::zero())
		{
			m_points.push_back(ss.position);
			m_colors.push_back(ss.irradiance / PI);
		}
		
	}

	return true;
}

void SceneSampleSet::addSample(SceneSample sample)
{
	m_samples.push_back(sample);
}

void SceneSampleSet::save()
{
	std::string savePath = "../data-files/Scenes/" + m_sceneName + "/SampleSets/" + m_sampleSetName + "/SamplePositions.txt";
	std::fstream saveFile;
	saveFile.open(savePath.c_str(), std::fstream::out);

	for (const SceneSample& ss : m_samples)
	{
		saveFile << ss.toString() << std::endl;
	}
	saveFile.close();
}

std::string generateUniqueName(int NumberOfProbes, G3D::Vector3 baseProbePosition)
{
	std::string baseLogPath = "C:/temp/logs/";
	std::stringstream ss;

#if WEIGHTS_CLOSEST
	int fileCount;
	Array<String> result;
	G3D::FileSystem::getFiles(baseLogPath.c_str(), result);
	fileCount = result.size();

	ss << "C:/temp/logs/";
	ss << "SS_";
	ss << NumberOfProbes << "probes_";
	ss << baseProbePosition.toString().c_str();
	ss << "bpp.txt";


#endif

	return ss.str();
}

//float SceneSampleSet::generateTriplet(int row, int col, const ProbeInterpolationRecord* iRec)
//{
//	float element = 0;
//
//	for (int coeff = 0; coeff < NumberOfCoefficientsPerProbe; ++coeff)
//	{
//		float gradient = probeStructure->getProbe(probeIndex)->coeffGradients[coeff][color][axis];
//
//		std::pair<int, int> lm = SH::kToLM(coeff);
//		float weight = iRec.weights[i];
//		float phong = phongCoeffs(lm.first, 1.0f);
//		float sh = SH::SHxyz_yup(lm.first, lm.second, SampleNormal);
//		element += gradient * phong * sh;
//	}
//	element *= iRec.weights[i];
//}

float SceneSampleSet::distanceToProbe(const G3D::Vector3& position, int p)
{
	Probe* probe = probeStructure->getProbe(p);
	return (position - probe->getPosition()).length();
}

float SceneSampleSet::SumOfInverseSquaredProbeDistances(const G3D::Vector3& position)
{
	float val = 0;
	for (int i = 0; i < probeStructure->probeCount(); ++i)
	{
		Probe* p = probeStructure->getProbe(i);
		val += 1.f / (powf((position - p->getPosition()).length(), 2));
	}
	return val;
}

float SceneSampleSet::InverseSumOf1OverSquaredProbeDistances(const G3D::Vector3& position)
{
	return 1.f / SumOfInverseSquaredProbeDistances(position);
}

float SceneSampleSet::dInverseSquaredSumdProbeN(const G3D::Vector3& position, const G3D::Vector3& normal, int m, int n, int axis, int color)
{
	return -1 * powf(InverseSumOf1OverSquaredProbeDistances(position), 2) * dInverseDistanceSquaredMdProbeN(position, normal, n, n, axis, color);
	//if (m != n)
	//{
	//	return 0;
	//}
	//Probe* pm = probeStructure->getProbe(m);
	//Probe* pn = probeStructure->getProbe(n);
	//const G3D::Vector3& posm = pm->getPosition();
	//const G3D::Vector3& posn = pn->getPosition();

	//const G3D::Vector3 posMinusPosProbeM = position - posm;
	//return -2 * powf(C(position), 2) * powf(posMinusPosProbeM.length(), -0.75) * (position[axis] - posn[axis]);
}

float SceneSampleSet::D(const G3D::Vector3& position, const G3D::Vector3& normal, int m, int n, int axis, int color)
{
	Probe* p = probeStructure->getProbe(m);
	return powf((position - p->getPosition()).length(), -2);

}

float SceneSampleSet::dInverseDistanceSquaredMdProbeN(const G3D::Vector3& position, const G3D::Vector3& normal, int m, int n, int axis, int color)
{
	if (m != n)
	{
		return 0;
	}

	Probe* pm = probeStructure->getProbe(m);
	Probe* pn = probeStructure->getProbe(n);
	const G3D::Vector3& posm = pm->getPosition();
	const G3D::Vector3& posn = pn->getPosition();

	const G3D::Vector3 posMinusPosProbeM = position - posm;


	return 2 * powf(posMinusPosProbeM.length(), -4) * (position[axis] - posm[axis]);

}

float SceneSampleSet::dWeightMdProbeN(const G3D::Vector3& position, const G3D::Vector3& normal, int m, int n, int axis, int color)
{
	return dInverseSquaredSumdProbeN(position, normal, m, n, axis, color) * powf(distanceToProbe(position, m), -2) + InverseSumOf1OverSquaredProbeDistances(position) * dInverseDistanceSquaredMdProbeN(position, normal, m, n, axis, color);
}

float SceneSampleSet::dRdX(const G3D::Vector3& position, const G3D::Vector3& normal, int NumberOfCoeffs, int m, int axis, int color)
{
#if defined(USE_CUSTOM_FUNCTION)
	return 1;
#else
	float value = 0;
#ifdef ONE_ROW_PER_SH_BAND
	float gradient = probeStructure->getProbe(m)->coeffGradients[NumberOfCoeffs][color][axis];
	std::pair<int, int> lm = SH::kToLM(NumberOfCoeffs);
	float sh = SH::SHxyz_yup(lm.first, lm.second, normal);
	float phong = phongCoeffs(lm.first, 1.0f);
	value += gradient /** phong * sh*/;

#else
	for (int coeff = 0; coeff < NumberOfCoeffs; ++coeff)
	{
		float gradient = probeStructure->getProbe(m)->coeffGradients[coeff][color][axis];
		
		std::pair<int, int> lm = SH::kToLM(coeff);
		float phong = phongCoeffs(lm.first, 1.0f);
		float sh = SH::SHxyz_yup(lm.first, lm.second, normal);
		value += gradient * phong * sh;
	}

#endif // ONE_ROW_PER_SH_BAND
	return value;
#endif
}

float SceneSampleSet::R(const G3D::Vector3& position,  const G3D::Vector3& normal, int NumberOfCoeffs,int m, int axis, int color)
{
#if defined(USE_CUSTOM_FUNCTION)
	return probeStructure->getProbe(m)->getPosition()[color];
#else
	float value = 0;
#ifdef ONE_ROW_PER_SH_BAND
	float coeffval = probeStructure->getProbe(m)->coeffs[NumberOfCoeffs][color];

	std::pair<int, int> lm = SH::kToLM(NumberOfCoeffs);
	float phong = phongCoeffs(lm.first, 1.0f);
	float sh = SH::SHxyz_yup(lm.first, lm.second, normal);
	value += coeffval /** phong * sh*/;
#else
	for (int coeff = 0; coeff < NumberOfCoeffs; ++coeff)
	{
		float coeffval = probeStructure->getProbe(m)->coeffs[coeff][color];

		std::pair<int, int> lm = SH::kToLM(coeff);
		float phong = phongCoeffs(lm.first, 1.0f);
		float sh = SH::SHxyz_yup(lm.first, lm.second, normal);
		value += coeffval * phong * sh;
	}
#endif // ONE_ROW_PER_SH_BAND
	return value;
#endif
}

float SceneSampleSet::w(const G3D::Vector3& position, const G3D::Vector3& normal, int m, int n, int axis, int color)
{
	return 0;
}


float SceneSampleSet::B(const G3D::Vector3& position,  const G3D::Vector3& normal, int NumberOfCoeffs,int m, int n, int axis, int color)
{
	if (m != n)
	{
		return 0;
	}
#if defined(USE_CUSTOM_FUNCTION)
	if (axis != color)
	{
		return 0;
	}
#endif
	return dRdX(position, normal, NumberOfCoeffs, m, axis, color);
}

void SceneSampleSet::generateTriplets(int NumberOfSamples,
                                      int NumberOfCoeffs,
                                      String outputPath, 
                                      std::vector<Eigen::Triplet<float>>* eigenTriplets, 
                                      bool optimizeForCoeffs)
{
    std::fstream outputFile;
    bool output = !(outputPath.empty());
    if (output)
    {
        outputFile.open(outputPath.c_str(), std::fstream::out);
    }
	
	int NumberOfProbes = probeStructure->probeCount();
	int NumberOfElements = NumberOfProbes * NumberOfCoeffs;

	//	// Scene (probe structure) Parameters
	int counter = 0;
	int row = 0;
	int col = 0;
	// Iterate over all wanted samples
	for (int sampleNumber = 0; sampleNumber < NumberOfSamples; ++sampleNumber)
	{
		//		bool increment = true;
		const SceneSample& ss = m_samples[sampleNumber];
		const G3D::Vector3& SamplePosition = ss.position;
		const G3D::Vector3& SampleNormal = invertNormal ? -ss.normal : ss.normal;

		ProbeInterpolationRecord iRec = probeStructure->getInterpolationProbeIndicesAndWeights(SamplePosition);
		//float weightDenum = C(SamplePosition);
		float inverseSumOf1OverSquaredProbeDistances = InverseSumOf1OverSquaredProbeDistances(SamplePosition);

#ifdef ONE_ROW_PER_SH_BAND
		for (int coeff = 0; coeff < NumberOfCoeffs; ++coeff)
		{
#endif

			for (int color = COLOR_RED; color <= COLOR_BLUE; ++color)
			{
				col = 0;



				for (int pn = 0; pn < NumberOfProbes; ++pn)
				{
					float dRGB_color_axis = 0;

					for (int axis = 0; axis <= AXIS_Z; ++axis)
					{
						int probeIndex = iRec.probeIndices[pn];
						float weight = iRec.weights[pn];

						Probe* n = probeStructure->getProbe(probeIndex);


						for (int pm = 0; pm < NumberOfProbes; ++pm)
						{
							Probe* m = probeStructure->getProbe(pm);

							float Wval = dWeightMdProbeN(SamplePosition, SampleNormal, pm, pn, axis, color);
#ifdef ONE_ROW_PER_SH_BAND
							float Rval = R(SamplePosition, SampleNormal, coeff, pm, axis, color);
							float Bval = B(SamplePosition, SampleNormal, coeff, pm, pn, axis, color);
#else
							float Rval = R(SamplePosition, SampleNormal, NumberOfCoeffs, pm, axis, color);
							float Bval = B(SamplePosition, SampleNormal, NumberOfCoeffs, pm, pn, axis, color);
#endif

							float computedWeights = inverseSumOf1OverSquaredProbeDistances * powf(distanceToProbe(SamplePosition, pm), -2);

							dRGB_color_axis += Wval * Rval + computedWeights * Bval;
						}

						if (eigenTriplets)
						{
							eigenTriplets->push_back(Eigen::Triplet<float>(row, col, dRGB_color_axis));
						}
						col++;
					}
				}
				row++;
			}
#ifdef ONE_ROW_PER_SH_BAND
		}
#endif
	}
}


void SceneSampleSet::generateRGBValuesFromProbes(int NumberOfSamples, int NumberOfCoeffs)
{
	m_colors.clear();

    std::fstream file = openFile(ESSFile::Values, false);
    for (int i = 0; i < NumberOfSamples; ++i)
    {
        const SceneSample& ss = m_samples[i];

        Vector3 rgb = probeStructure->reconstructSH(ss.position, ss.normal, NumberOfCoeffs);

        file << rgb.x << " " << rgb.y << " " << rgb.z << "\n";
		m_colors.push_back(Color3(rgb));
	}

    file.close();
}

void SceneSampleSet::generateInterpolatedCoefficientsFromProbes(int NumberOfSamples, int NumberOfCoeffs, String savePath, Eigen::VectorXd* eigenVector)
{
	//	// Scene (probe structure) Parameters
	std::fstream samplesRGBFile;
	std::fstream logFile;

	samplesRGBFile.open(savePath.c_str(), std::ios::out);
	samplesRGBFile.precision(20);

	for (int i = 0; i < NumberOfSamples; ++i)
	{
		const SceneSample& ss = m_samples[i];

		TProbeCoefficients interpolatedCoeffs = probeStructure->interpolatedCoefficients(ss.position, ss.normal, NumberOfCoeffs);

		dumpToFile(samplesRGBFile, interpolatedCoeffs);

		if (eigenVector)
		{
			for (int j = 0; j < NumberOfCoeffs; ++j)
			{
				(*eigenVector)(i * (3 * NumberOfCoeffs) + j * 3 + 0) = interpolatedCoeffs[j][0];
				(*eigenVector)(i * (3 * NumberOfCoeffs) + j * 3 + 1) = interpolatedCoeffs[j][1];
				(*eigenVector)(i * (3 * NumberOfCoeffs) + j * 3 + 2) = interpolatedCoeffs[j][2];
			}
		}
	}
}

void SceneSampleSet::generateRGBValuesFromSamples(int NumberOfSamples, String savePath, Eigen::VectorXd* eigenVector)
{
    std::fstream samplesRGBFile(savePath.c_str(), std::ios::out);;
    samplesRGBFile.precision(20);

    for (int i = 0; i < NumberOfSamples; ++i)
    {
        const SceneSample& ss = m_samples[i];
        Vector3 rgb = Vector3(ss.irradiance);
        samplesRGBFile << rgb.x << std::endl;
        samplesRGBFile << rgb.y << std::endl;
        samplesRGBFile << rgb.z << std::endl;

        if (eigenVector)
        {
            (*eigenVector)(i * 3 + 0) = rgb.x;
            (*eigenVector)(i * 3 + 1) = rgb.y;
            (*eigenVector)(i * 3 + 2) = rgb.z;
        }
    }

    samplesRGBFile.close();
}

void SceneSampleSet::generateRGBValuesFromProbes(int NumberOfSamples, int NumberOfCoeffs, String savePath, Eigen::VectorXd* eigenVector = 0)
{
	int NumberOfProbes = probeStructure->probeCount();
	int NumberOfElements = NumberOfProbes * NumberOfCoeffs;

	//	// Scene (probe structure) Parameters
    std::fstream samplesRGBFile;
    bool output = !(savePath.empty());
    if (output)
    {
        samplesRGBFile.open(savePath.c_str(), std::ios::out);;
        samplesRGBFile.precision(20);
    }

    for (int i = 0; i < NumberOfSamples; ++i)
	{
        const SceneSample& ss = m_samples[i];

#if defined(USE_CUSTOM_FUNCTION)
		ProbeInterpolationRecord record = probeStructure->getInterpolationProbeIndicesAndWeights(ss.position);

		Vector3 rgb;
		for (int i = 0; i < record.weights.size(); ++i)
		{
			rgb += record.weights[i] * probeStructure->getProbe(i)->getPosition();
		}

#else
		
		Vector3 rgb = probeStructure->reconstructSH(ss.position, ss.normal, NumberOfCoeffs);
#endif


        if (output)
        {
            samplesRGBFile << rgb.x << std::endl;
            samplesRGBFile << rgb.y << std::endl;
            samplesRGBFile << rgb.z << std::endl;
        }

		if (eigenVector)
		{
			(*eigenVector)(i * 3 + 0) = rgb.x;
			(*eigenVector)(i * 3 + 1) = rgb.y;
			(*eigenVector)(i * 3 + 2) = rgb.z;
		}
	}
}

void SceneSampleSet::createbVector(Eigen::VectorXd* bVector, const Eigen::VectorXd* rgbColumn, String& optimizationFolderPath)
{
	// this should not pull from a file but rather straight up use the values in the scenesamples
	std::fstream refFile;
	refFile.open((optimizationFolderPath + "/ref_values.txt").c_str(), std::ios::in);
	std::string line;
	int counter = 0;
	while (std::getline(refFile, line))
	{
		float value = std::stof(line.c_str());
		float value2 = (*rgbColumn)(counter);
		float toWrite = value - value2;
		(*bVector)(counter) = toWrite;
		counter++;
	}
	refFile.close();
}

std::fstream SceneSampleSet::openFile(ESSFile type, bool reading)
{
    std::string val;
    if (type == ESSFile::Samples)
    {
        val = "SamplePositions.txt";
    }
    else
    {
        val = "IrradianceResults2.txt";
    }

    int operation;

    if (reading)
    {
        operation = std::fstream::in;
    }
    else
    {
        operation = std::fstream::out;
    }

    std::string path = "../data-files/Scenes/" + m_sceneName + "/SampleSets/" + m_sampleSetName + "/" + val;
    return std::fstream(path.c_str(), operation);
}

void SceneSampleSet::clearValues()
{
	std::fstream irradianceFile;
	
	std::string irradiancePath = "../data-files/Scenes/" + m_sceneName + "/SampleSets/" + m_sampleSetName + "/IrradianceResults2.txt";

	irradianceFile.open(irradiancePath.c_str(), std::fstream::out);
	irradianceFile.close();
	m_colors.clear();
}

void SceneSampleSet::clearPositions()
{
	std::fstream samplesFile;
	std::string samplesPath = "../data-files/Scenes/" + m_sceneName + "/SampleSets/" + m_sampleSetName + "/SamplePositions.txt";


	samplesFile.open(samplesPath.c_str(), std::fstream::out);
	samplesFile.close();
	m_points.clear();
}

WeightMatrixType SceneSampleSet::generateWeightsMatrix(int NumberOfSamples, int NumberOfCoeffs)
{
    int NumElementsPerRow = probeStructure->probeCount() * 3;

    std::vector<Eigen::Triplet<float>>* eigenTriplets = new std::vector<Eigen::Triplet<float>>;
    eigenTriplets->reserve(NumberOfSamples * 3 * NumElementsPerRow);

    generateTriplets(NumberOfSamples, NumberOfCoeffs, "", eigenTriplets, false);

    WeightMatrixType A(NumberOfSamples * 3, NumElementsPerRow);
    A.setFromTriplets(eigenTriplets->begin(), eigenTriplets->end());

    return A;
}

Eigen::VectorXd SceneSampleSet::generatebVector(int NumberOfSamples, int NumberOfCoeffs, String optimizationFolderPath)
{
	Eigen::VectorXd* rgbColumn = new Eigen::VectorXd(NumberOfSamples * 3);
	generateRGBValuesFromProbes(NumberOfSamples, NumberOfCoeffs, "", rgbColumn);

	Eigen::VectorXd bVector(NumberOfSamples * 3);
	createbVector(&bVector, rgbColumn, optimizationFolderPath);

	return bVector;
}

void SceneSampleSet::outputWeightsMatrixToFile(int NumberOfSamples, int NumberOfCoeffs, String optimizationFolderPath)
{
    std::fstream outputFile((optimizationFolderPath + "/A.txt").c_str(), std::fstream::out);

    WeightMatrixType A = generateWeightsMatrix(NumberOfSamples, NumberOfCoeffs);
    outputFile << A;

    outputFile.close();
}

void SceneSampleSet::outputBVectorToFile(int NumberOfSamples, int NumberOfCoeffs, String optimizationFolderPath)
{
    std::fstream outputFile((optimizationFolderPath + "/b.txt").c_str(), std::fstream::out);

    Eigen::VectorXd* rgbColumn = new Eigen::VectorXd(NumberOfSamples * 3);
    generateRGBValuesFromProbes(NumberOfSamples, NumberOfCoeffs, "", rgbColumn);

    Eigen::VectorXd bVector(NumberOfSamples * 3);
    createbVector(&bVector, rgbColumn, optimizationFolderPath);

	outputFile << bVector;
    outputFile.close();
}

std::vector<float> SceneSampleSet::tryOptimizationPass(int NumberOfSamples, 
													   int NumberOfCoeffs, 
													   bool optimizeForMitsubaSamples, 
													   String optimizationFolderPath)
{
	G3D::StopWatch sw("InnerOptimization");
	sw.setEnabled(true);

	// todo: this is dependant on the interpolation method...
	int NumElementsPerRow = probeStructure->probeCount() * 3;

	std::vector<Eigen::Triplet<float>>* eigenTriplets = new std::vector<Eigen::Triplet<float>>;
	eigenTriplets->reserve(NumberOfSamples * 3 * NumElementsPerRow);


	Eigen::VectorXd* rgbColumn = new Eigen::VectorXd(NumberOfSamples * 3 * NumberOfCoeffs);

	generateTriplets(NumberOfSamples, NumberOfCoeffs, optimizationFolderPath + "/triplets.txt", eigenTriplets, optimizeForMitsubaSamples);
	sw.after("Generated Triplets");
	//generateRGBValuesFromProbes(NumberOfSamples, NumberOfCoeffs, optimizationFolderPath + "/values.txt", rgbColumn);
	generateInterpolatedCoefficientsFromProbes(NumberOfSamples, NumberOfCoeffs, optimizationFolderPath + "/values.txt", rgbColumn);
	
	Eigen::VectorXd bVector(NumberOfSamples * 3 * NumberOfCoeffs);
	createbVector(&bVector, rgbColumn, optimizationFolderPath);

	WeightMatrixType A(NumberOfSamples * 3 * NumberOfCoeffs, NumElementsPerRow);
	A.setFromTriplets(eigenTriplets->begin(), eigenTriplets->end());
	sw.after("Set matrices");

	int NumberOfProbes = probeStructure->probeCount();

	Eigen::VectorXd optimizationResult(NumberOfProbes * 3);
	bool success = probeOptimizationPass(A, bVector, &optimizationResult);
	sw.after("Optimization finished");

	std::fstream logFile;
	logFile.open((optimizationFolderPath + "/log.txt").c_str(), std::ios::app);
	logFile << "After optimization \n";
	logFile << "A\n" << A << "\n";
	logFile << "b\n" << bVector << "\n";
    logFile << "x\n" << optimizationResult << "\n";
    logFile << "A*x\n" << A*optimizationResult << "\n";
	logFile.close();

	std::vector<float> toReturn;
	if (success)
	{
		for (int i = 0; i < NumberOfProbes * 3; ++i)
		{
			toReturn.push_back((float)optimizationResult(i));
		}
		
		std::fstream outputFile;
		outputFile.open((optimizationFolderPath + "/dp.txt").c_str(), std::ios::out);
		outputFile.precision(20);
		outputFile << optimizationResult;
		outputFile.close();

		std::fstream dpLogFile;
		dpLogFile.open((optimizationFolderPath + "/dplog.txt").c_str(), std::fstream::out | std::fstream::app | std::fstream::in);

		if (!dpLogFile.is_open())
		{
			// create file because it does not exist
			dpLogFile.open((optimizationFolderPath + "/dplog.txt").c_str(),  std::fstream::in | std::fstream::out | std::fstream::trunc);
		}

        dpLogFile << optimizationResult << "\n\n";
	}
	sw.after("Logged results");


	delete eigenTriplets;
	return toReturn;
}


bool SceneSampleSet::probeOptimizationPass(WeightMatrixType& A, Eigen::VectorXd& b, Eigen::VectorXd* result)
{
	WeightMatrixType At = A.transpose();
	WeightMatrixType AtA = At * A;
	Eigen::VectorXd Atb = At * b;


	Eigen::ConjugateGradient<WeightMatrixType> solver;
	debugPrintf("Compute step...\n");
	solver.setTolerance(1e-4);
	solver.setMaxIterations(1e12);
	solver.compute(AtA);

	if (solver.info() != Eigen::Success)
	{
		debugPrintf("Compute step failed!");
		return false;
	}
	std::fstream logFile;
	logFile.open("C:/temp/log.txt", std::ios::out);
	//logFile << A << "\n";
	//logFile << At << "\n";
	//logFile << AtA << "\n";
	logFile << "b\n" << b << "\n\n";
	logFile << "Atb\n" << Atb << "\n\n";
	debugPrintf("Solve step...\n");
	Eigen::VectorXd tempResult = solver.solve(Atb);
	//
	logFile << "Result\n" << AtA * tempResult << "\n\n";
	logFile << "Result2\n" << A * tempResult;
	logFile.close();
	//
	for (int i = 0; i < tempResult.size(); ++i)
	{
		//debugPrintf("%f\n", tempResult(i));
		(*result)(i) = tempResult(i);
	}
	if (solver.info() != Eigen::Success)
	{
		debugPrintf("Solve step FAILED!\n");
		return false;
	}

	return true;
}