#include <iostream>
#include <fstream>
//#include <utility>
#include <string>
#include <sstream>

#include "SceneSampleSet.h"
#include "SceneSample.h"
#include "ProbeStructure.h"
#include "SH.h"


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

#define PI 3.141592653589793238462643383279


bool invertNormal = true;
void probeOptimizationPass(std::vector<int> matrixRows, std::vector<int> matrixColumns, std::vector<float> matrixElements, std::vector<float> rgb);

SceneSampleSet::SceneSampleSet(std::string sceneName, std::string sampleSetName, float scale)
{
	m_sceneName = sceneName;
	m_sampleSetName = sampleSetName;
	m_scale = scale;
	load();
}
void SceneSampleSet::load()
{
	m_samples.clear();
	std::fstream samplesFile, irradianceFile;
	std::string samplesPath = "C:/libraries/g3d/samples/aTest/data-files/Scenes/" + m_sceneName + "/SampleSets/" + m_sampleSetName + "/SamplePositions.txt";

	bool loadIrradianceResults2 = true;
	std::string irradiancePath;

	if (loadIrradianceResults2)
	{
		irradiancePath = "C:/libraries/g3d/samples/aTest/data-files/Scenes/" + m_sceneName + "/SampleSets/" + m_sampleSetName + "/IrradianceResults2.txt";
	}
	else
	{
		irradiancePath = "C:/libraries/g3d/samples/aTest/data-files/Scenes/" + m_sceneName + "/SampleSets/" + m_sampleSetName + "/IrradianceResults.txt";

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

			//if (loadIrradianceResults2)
			//{
			//	ss.irradiance = Color3(std::atof(splitLine[0].c_str()),
			//						   std::atof(splitLine[1].c_str()),
			//						   std::atof(splitLine[2].c_str()));
			//	//ss.irradiance *= 10;
			//}
			//else
			//{
				ss.irradiance = Color3(std::atof(splitLine[0].c_str()) / 255.f,
					std::atof(splitLine[1].c_str()) / 255.f,
					std::atof(splitLine[2].c_str()) / 255.f);
			//}

		}

		if (currentLine == 2)
		{ // Samples Positions
			ss.position = Vector3(std::atof(splitLine[0].c_str()) * m_scale,
								  std::atof(splitLine[1].c_str()) * m_scale,
								  std::atof(splitLine[2].c_str()) * m_scale);
		}

		if (currentLine == 4)
		{ // Normals
			ss.normal = Vector3(std::atof(splitLine[0].c_str()),
				std::atof(splitLine[1].c_str()),
				std::atof(splitLine[2].c_str()));
		}

		currentLine++;

		if (currentLine == NumberOfLinesPerSample)
		{
			m_samples.push_back(ss);
			currentLine = 0;
		}

	}

	samplesFile.close();
	irradianceFile.close();

}


void SceneSampleSet::addSample(SceneSample sample)
{
	m_samples.push_back(sample);
}
void SceneSampleSet::save()
{
	std::string savePath = "C:/libraries/G3D/samples/aTest/data-files/Scenes/" + m_sceneName + "/SampleSets/" + m_sampleSetName + "/SamplePositions.txt";
	std::fstream saveFile;
	saveFile.open(savePath.c_str(), std::fstream::out);

	for (const SceneSample& ss : m_samples)
	{
		saveFile << ss.toString() << std::endl;
	}
	saveFile.close();
}

double phongCoeffs(int l, double r = 1.0f)
{
	if (l == 0)
	{
		return PI;
	}
	else if (l == 1)
	{
		return PI * (1.0 + r) / (2.0 + r);
	}
	else if (l == 2)
	{
		return PI * r / (3.0 + r);
	}
	return 0;
}

double NDotOmegaCoeff(int l)
{
	if (l == 1)
	{
		return 2 * PI / 3;
	}
	else if (l % 2 == 0)
	{
		int a = pow(-1, l / 2.f - 1);
		int b = SH::factorial(l);
		return (2 * PI * a * b) / ((l + 2) * (l + 1) * pow(2, l) * SH::factorial(pow(l / 2.f, 2)));
	}
	else
	{
		return 0;
	}
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


void SceneSampleSet::generateTriplets(int NumberOfSamples, std::vector<Eigen::Triplet<float>>* eigenTriplets = 0)
{
	std::fstream outputFile;
	outputFile.open("C:/temp/triplets2.txt", std::fstream::out);

	int NumberOfCoefficientsPerProbe = 9;
	int NumberOfProbes = probeStructure->probeCount();
	int NumberOfElements = NumberOfProbes * NumberOfCoefficientsPerProbe;

	//	// Scene (probe structure) Parameters
	const G3D::Vector3 firstProbePosition = G3D::Vector3(probeStructure->m_firstProbePosition[0], probeStructure->m_firstProbePosition[1], probeStructure->m_firstProbePosition[2]);
	const G3D::Vector3 dimensions = G3D::Vector3(probeStructure->m_dimensions[0], probeStructure->m_dimensions[1], probeStructure->m_dimensions[2]);
	float step = probeStructure->m_step;
	int counter = 0;
	int elementCounter = 0;
	// Iterate over all wanted samples
	for (int sampleNumber = 0; sampleNumber < NumberOfSamples; ++sampleNumber)
	{
		//		bool increment = true;
		const SceneSample& ss = m_samples[sampleNumber];
		//std::cout << ss.toString() << std::endl;
		const G3D::Vector3& SamplePosition = ss.position;
		const G3D::Vector3& SampleNormal = invertNormal? -ss.normal : ss.normal;


		ProbeInterpolationRecord iRec = probeStructure->getInterpolationProbeIndicesAndWeights(SamplePosition);

		// For all interpolated probes
		for (int i = 0; i < iRec.weights.size(); ++i)
		{
			//const G3D::Vector3& probeCoords = coords[i];

#if GENERATE_POSITION_FITTING_MATRIX

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
						float phong = phongCoeffs(lm.first);
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

#endif
#if GENERATE_COEFFICIENTS_FITTING_MATRIX
			//			// If the probe is outside the structure, we don't consider it
			//
			//
			//			// Find the column at which we start entering factors
			int startIndex = getProbeStartIndex(probeCoords, firstProbePosition, dimensions, NumberOfCoefficientsPerProbe);

			if ((startIndex < 0) || (startIndex >= NumberOfElements))
			{
				continue;
			}

#if LIMIT_CELL_OCCUPANCY
			if (i == 0)
			{
				//std::cout << startIndex / 9 << std::endl;
				cellCounts[startIndex / 9.f] += 1;

				if (cellCounts[startIndex / 9.f] >= MAX_CELL_OCCUPANCY)
				{
					increment = false;
					break;
				}
			}
#endif

			// For all SH bands
			for (int coeff = 0; coeff < NumberOfCoefficientsPerProbe; ++coeff)
			{
				int* lm = kToLM(coeff);

				int col = startIndex + coeff;

				// We multiply the inter weight by the geometric term and the SH function value in this direction
				double factors = interpolationWeights[i] *
					//NDotOmegaCoeff(lm(0)) *
					phongCoeffs(lm[0]) *
					SHxyz_yup(lm[0], lm[1], SampleNormal);

				//triplets.push_back(Eigen::Triplet<double>(counter, col, factors));
				outputFile << counter << " " << col << " " << factors << std::endl;

			}
#endif
		}
		//
		//#if LIMIT_CELL_OCCUPANCY
		//		if (increment)
		//#endif
		{
			//counter++;
		}

	}
	outputFile << counter << " " << probeStructure->probeCount()*3-1 << " " << 0 << std::endl;

}

void SceneSampleSet::generateRGBValuesFromProbes(int NumberOfSamples, bool ref = false, Eigen::VectorXd* eigenVector = 0)
{

	int NumberOfCoefficientsPerProbe = 9;
	int NumberOfProbes = probeStructure->probeCount();
	int NumberOfElements = NumberOfProbes * NumberOfCoefficientsPerProbe;

	//	// Scene (probe structure) Parameters
	const G3D::Vector3 firstProbePosition = G3D::Vector3(probeStructure->m_firstProbePosition[0], probeStructure->m_firstProbePosition[1], probeStructure->m_firstProbePosition[2]);
	const G3D::Vector3 dimensions = G3D::Vector3(probeStructure->m_dimensions[0], probeStructure->m_dimensions[1], probeStructure->m_dimensions[2]);
	float step = probeStructure->m_step;


	std::fstream samplesRGB;
	std::fstream logFile;

	if (ref)
	{
		samplesRGB.open("C:/temp/samplesRGB_ref.txt", std::ios::out);
	}
	else
	{
		samplesRGB.open("C:/temp/samplesRGB.txt", std::ios::out);
	}
	samplesRGB.precision(20);

	int counter = 0;
	if (this->outputToLog)
	{
		std::string logFileName = generateUniqueName(NumberOfProbes, firstProbePosition);
		logFile.open(logFileName.c_str(), std::ios::out);
	}

	for (SceneSample sample : m_samples)
	{
		const G3D::Vector3& SamplePosition = sample.position;
		const G3D::Vector3& SampleNormal = invertNormal? -sample.normal : sample.normal;

		ProbeInterpolationRecord iRec = probeStructure->getInterpolationProbeIndicesAndWeights(SamplePosition);

		if (this->outputToLog)
		{
			logFile << "Position = " << SamplePosition.toString().c_str() << std::endl;
			logFile << "Normal = " << SampleNormal.toString().c_str() << std::endl << std::endl;
		}

		Vector3 rgb = Vector3(0, 0, 0);
		for (int i = 0; i < iRec.weights.size(); ++i)
		{
			//int probeIndex = getProbeStartIndex(coords[i], firstProbePosition, dimensions, 1);
			int probeIndex = iRec.probeIndices[i];
			float weight = iRec.weights[i];

			if (this->outputToLog)
			{
				logFile << "Probe " << i << std::endl;
				//logFile << "position =  " << coords[i].toString().c_str() << std::endl;
				logFile << "weight =  " << weight << std::endl;
			}

			// For all SH bands
			for (int coeff = 0; coeff < NumberOfCoefficientsPerProbe; ++coeff)
			{
				std::pair<int, int> lm = SH::kToLM(coeff);

				// We multiply the inter weight by the geometric term and the SH function value in this direction
				float phong = phongCoeffs(lm.first);
				float sh = SH::SHxyz_yup(lm.first, lm.second, SampleNormal);
				double factors = weight *
					//NDotOmegaCoeff(lm(0)) *
					phong *
					sh;
				Vector3& probeCoeffs = probeStructure->getProbe(probeIndex)->coeffs[coeff];
				Vector3 accumulation = factors * probeCoeffs;
				rgb += accumulation;

				if (this->outputToLog)
				{
					logFile << "SH band: ( " << lm.first << ", " << lm.second << " )" << std::endl;
					logFile << "phong = " << phong << std::endl;
					logFile << "shFunctionEvaluation = " << sh << std::endl;
					logFile << "weight * phong * sh = " << factors << std::endl;
					logFile << "probeCoeffs (R,G,B) = " << probeCoeffs.toString().c_str() << std::endl;
					logFile << "weight * phong * sh * probeCoeffs = " << accumulation.toString().c_str() << std::endl;
					logFile << "accumulatedRGB(" << coeff << ") = " << rgb.toString().c_str() << std::endl << std::endl;
				}

			}
		}
		rgb /= 3.141592654;
		samplesRGB << rgb.x << std::endl;
		samplesRGB << rgb.y << std::endl;
		samplesRGB << rgb.z << std::endl;

		if (eigenVector)
		{
			(*eigenVector)(counter * 3 + 0) = rgb.x;
			(*eigenVector)(counter * 3 + 1) = rgb.y;
			(*eigenVector)(counter * 3 + 2) = rgb.z;
		}

		if (++counter == NumberOfSamples)
		{
			break;
		}
	}
}

void SceneSampleSet::createbVector(Eigen::VectorXd* bVector, const Eigen::VectorXd* rgbColumn)
{
	// this should not pull from a file but rather straight up use the values in the scenesamples
	std::fstream refFile;
	refFile.open("C:/temp/samplesRGB_ref.txt", std::ios::in);
	std::string line;
	int counter = 0;
	while (std::getline(refFile, line))
	{
		float value = std::atof(line.c_str());
		(*bVector)(counter) = value - (*rgbColumn)(counter);
		counter++;
	}
	refFile.close();
}

std::vector<float> SceneSampleSet::tryOptimizationPass(int NumberOfSamples, bool ref)
{
	
	// todo: this is dependant on the interpolation method...
	int NumElementsPerRow = probeStructure->probeCount() * 3;

	std::vector<Eigen::Triplet<float>>* eigenTriplets = new std::vector<Eigen::Triplet<float>>;
	eigenTriplets->reserve(NumberOfSamples * 3 * NumElementsPerRow);


	Eigen::VectorXd* rgbColumn = new Eigen::VectorXd(NumberOfSamples * 3);

	generateTriplets(NumberOfSamples, eigenTriplets);
	generateRGBValuesFromProbes(NumberOfSamples, ref, rgbColumn);

	Eigen::VectorXd bVector(NumberOfSamples * 3);
	createbVector(&bVector, rgbColumn);

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
			toReturn.push_back(optimizationResult(i));
		}
		
		std::fstream outputFile;
		outputFile.open("C:/temp/dp.txt", std::ios::out);
		outputFile.precision(20);
		outputFile << optimizationResult;
		outputFile.close();

		std::fstream dpLogFile;
		dpLogFile.open("C:/temp/CurrentOptimization/dplog.txt", std::fstream::out | std::fstream::app | std::fstream::in);

		if (!dpLogFile)
		{
			// create file because it does not exist
			dpLogFile.open("C:/temp/CurrentOptimization/dplog.txt",  std::fstream::in | std::fstream::out | std::fstream::trunc);
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


	//Eigen::LeastSquaresConjugateGradient<WeightMatrixType> solver;
	//debugPrintf("Compute step...\n");
	//solver.setTolerance(1e-5);
	//solver.setMaxIterations(1e5);
	//solver.compute(AtA);

	//if (solver.info() != Eigen::Success)
	//{
	//	debugPrintf("Compute step failed!");
	//	return false;
	//}
	//std::fstream logFile;
	//logFile.open("C:/temp/log.txt", std::ios::out);
	////logFile << A << "\n";
	////logFile << At << "\n";
	////logFile << AtA << "\n";
	//logFile << "b\n" << b << "\n\n";
	//logFile << "Atb\n" << Atb << "\n\n";
	//debugPrintf("Solve step...");
	//Eigen::VectorXd tempResult = solver.solve(Atb);
	////
	//logFile << "Result\n" << AtA * tempResult << "\n\n";
	//logFile << "Result2\n" << A * tempResult;
	//logFile.close();
	////
	//for (int i = 0; i < tempResult.size(); ++i)
	//{
	//	//debugPrintf("%f\n", tempResult(i));
	//	(*result)(i) = tempResult(i);
	//}
	//if (solver.info() != Eigen::Success)
	//{
	//	debugPrintf("Solve step failed!");
	//	return false;
	//}

	return true;
}