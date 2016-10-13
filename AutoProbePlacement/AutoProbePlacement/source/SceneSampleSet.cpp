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


bool invertNormal = true;
void probeOptimizationPass(std::vector<int> matrixRows, std::vector<int> matrixColumns, std::vector<float> matrixElements, std::vector<float> rgb);


SceneSampleSet::SceneSampleSet(std::string sceneName, std::string sampleSetName, float scale, int numSamplesToLoad)
{
	m_sceneName = sceneName;
	m_sampleSetName = sampleSetName;
	m_scale = scale;
	load(numSamplesToLoad);
}
void SceneSampleSet::load(int maxSamples)
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
		return;
	}

	if (!irradianceFile.is_open())
	{
		debugPrintf("Couldn't open irradiance results file!\n");
		return;
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
				break;
			}

			Array<String> splitLine = stringSplit(String(irradianceLine.c_str()), ' ');

			float divider = 1.f;
			ss.irradiance = Color3(std::stof(splitLine[0].c_str()) / divider,
								   std::stof(splitLine[1].c_str()) / divider,
								   std::stof(splitLine[2].c_str()) / divider);
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
		m_points.push_back(ss.position);
		m_colors.push_back(ss.irradiance / PI);
	}
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

void SceneSampleSet::generateTriplets(int NumberOfSamples, String outputPath, std::vector<Eigen::Triplet<float>>* eigenTriplets = 0)
{
	std::fstream outputFile;
	outputFile.open(outputPath.c_str(), std::fstream::out);

	int NumberOfCoefficientsPerProbe = 9;
	int NumberOfProbes = probeStructure->probeCount();
	int NumberOfElements = NumberOfProbes * NumberOfCoefficientsPerProbe;

	//	// Scene (probe structure) Parameters
	int counter = 0;

	// Iterate over all wanted samples
	for (int sampleNumber = 0; sampleNumber < NumberOfSamples; ++sampleNumber)
	{
		//		bool increment = true;
		const SceneSample& ss = m_samples[sampleNumber];
		const G3D::Vector3& SamplePosition = ss.position;
		const G3D::Vector3& SampleNormal = invertNormal? -ss.normal : ss.normal;

		ProbeInterpolationRecord iRec = probeStructure->getInterpolationProbeIndicesAndWeights(SamplePosition);

		// For all interpolated probes
		for (int i = 0; i < iRec.weights.size(); ++i)
		{
			int probeIndex = iRec.probeIndices[i];
			int startIndex = probeIndex * 3;
			for (int color = COLOR_RED; color <= COLOR_BLUE; ++color)
			{
				counter = sampleNumber * 3 + color;
				for (int axis = AXIS_X; axis <= AXIS_Z; ++axis)
				{
					float element = 0;

					for (int coeff = 0; coeff < NumberOfCoefficientsPerProbe; ++coeff)
					{
						float gradient = probeStructure->getProbe(probeIndex)->coeffGradients[coeff][color][axis];

						std::pair<int, int> lm = SH::kToLM(coeff);
						float weight = iRec.weights[i]; 
						float phong = phongCoeffs(lm.first, 1.0f);
						float sh = SH::SHxyz_yup(lm.first, lm.second, SampleNormal);
						element += gradient * phong * sh;
					}
					element *= iRec.weights[i];
					int col = startIndex + axis;
					outputFile << counter << " " << col << " " << element << std::endl;
					if (eigenTriplets)
					{
						eigenTriplets->push_back(Eigen::Triplet<float>(counter, col, element));
					}
				}
			}
		}
	}
	outputFile << counter << " " << probeStructure->probeCount()*3-1 << " " << 0 << std::endl;
}

void SceneSampleSet::generateRGBValuesFromProbes(int NumberOfSamples)
{
    std::fstream file = openFile(ESSFile::Values, false);
    for (int i = 0; i < NumberOfSamples; ++i)
    {
        const SceneSample& ss = m_samples[i];

        Vector3 rgb = probeStructure->reconstructSH(ss.position, ss.normal);

        file << rgb.x << " " << rgb.y << " " << rgb.z << "\n";
    }

    file.close();
}

void SceneSampleSet::generateRGBValuesFromProbes(int NumberOfSamples, String savePath, Eigen::VectorXd* eigenVector = 0)
{
	int NumberOfCoefficientsPerProbe = 9;
	int NumberOfProbes = probeStructure->probeCount();
	int NumberOfElements = NumberOfProbes * NumberOfCoefficientsPerProbe;

	//	// Scene (probe structure) Parameters
	std::fstream samplesRGBFile;
	std::fstream logFile;

    samplesRGBFile.open(savePath.c_str(), std::ios::out);
    samplesRGBFile.precision(20);

    for (int i = 0; i < NumberOfSamples; ++i)
	{
        const SceneSample& ss = m_samples[i];

        Vector3 rgb = probeStructure->reconstructSH(ss.position, ss.normal);

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
}

void SceneSampleSet::createbVector(Eigen::VectorXd* bVector, const Eigen::VectorXd* rgbColumn, String& optimizationFolderPath)
{
	// this should not pull from a file but rather straight up use the values in the scenesamples
	std::fstream refFile;
	refFile.open((optimizationFolderPath + "/samplesRGB_ref.txt").c_str(), std::ios::in);
	std::string line;
	int counter = 0;
	while (std::getline(refFile, line))
	{
		float value = std::stof(line.c_str());
		(*bVector)(counter) = value - (*rgbColumn)(counter);
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




std::vector<float> SceneSampleSet::tryOptimizationPass(int NumberOfSamples, bool ref, String optimizationFolderPath)
{
	// todo: this is dependant on the interpolation method...
	int NumElementsPerRow = probeStructure->probeCount() * 3;

	std::vector<Eigen::Triplet<float>>* eigenTriplets = new std::vector<Eigen::Triplet<float>>;
	eigenTriplets->reserve(NumberOfSamples * 3 * NumElementsPerRow);


	Eigen::VectorXd* rgbColumn = new Eigen::VectorXd(NumberOfSamples * 3);

	generateTriplets(NumberOfSamples, optimizationFolderPath + "/triplets.txt", eigenTriplets);
	generateRGBValuesFromProbes(NumberOfSamples, optimizationFolderPath + "/samplesRGB.txt", rgbColumn);

	Eigen::VectorXd bVector(NumberOfSamples * 3);
	createbVector(&bVector, rgbColumn, optimizationFolderPath);

	WeightMatrixType A(NumberOfSamples * 3, NumElementsPerRow);
	A.setFromTriplets(eigenTriplets->begin(), eigenTriplets->end());

	int NumberOfProbes = probeStructure->probeCount();

	Eigen::VectorXd optimizationResult(NumberOfProbes * 3);
	bool success = probeOptimizationPass(A, bVector, &optimizationResult);

	//std::fstream logFile;
	//logFile.open("C:/temp/log.txt", std::ios::app);
	//logFile << "After optimization \n";
	//logFile << A << "\n";
	//logFile << bVector << "\n";
	//logFile << optimizationResult << "\n";
	//logFile.close();

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

		if (!dpLogFile)
		{
			// create file because it does not exist
			dpLogFile.open((optimizationFolderPath + "/dplog.txt").c_str(),  std::fstream::in | std::fstream::out | std::fstream::trunc);
		}
		dpLogFile << optimizationResult;
	}


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
	solver.setTolerance(1e-5);
	solver.setMaxIterations(1e5);
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
	debugPrintf("Solve step...");
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
		debugPrintf("Solve step failed!");
		return false;
	}

	return true;
}