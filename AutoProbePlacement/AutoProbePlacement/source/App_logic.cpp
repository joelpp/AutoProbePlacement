#include "App.h"
#include "SceneSampleSet.h"
#include "Helpers.h"

Vector3 App::getRandomPointInScene(){
	float x = Random::common().uniform(-5,5);
	float y = Random::common().uniform(0,5);
	float z = Random::common().uniform(-5,5);

	return Vector3(x,y,z);
}

Array<Vector3> App::getRandomPoint(int modelNumber, Vector3* P, Vector3* N, Vector3* barycentricWeights, int* s){
    Array<Vector3> toReturn = Array<Vector3>();
	//Random r(System::time());

	MeshShape* meshShape = m_scene->getMeshShape(modelNumber);

    meshShape->getRandomSurfacePoint(*P, *N, *s, *barycentricWeights, Random::common());

    //debugPrintf("triangle_vertex1: %s\n",meshShapes[modelNumber]->vertexArray()[*s%4].toString().c_str());
    //debugPrintf("triangle_vertex2: %s\n",meshShapes[modelNumber]->vertexArray()[(*s+1)%4].toString().c_str());
    //debugPrintf("triangle_vertex3: %s\n",meshShapes[modelNumber]->vertexArray()[(*s+2)%4].toString().c_str());

    int s1 = meshShape->indexArray()[*s];
    int s2 = meshShape->indexArray()[*s+1];
    int s3 = meshShape->indexArray()[*s+2];

    toReturn.append(Vector3(meshShape->vertexArray()[s1]));
    toReturn.append(Vector3(meshShape->vertexArray()[s2]));
    toReturn.append(Vector3(meshShape->vertexArray()[s3]));
    randomPoints.append(*P);
    return toReturn;
}

SceneSample App::generateSceneSample()
{
    Vector3* baryWeights = new Vector3();
    Vector3* P = new Vector3();
    Vector3* N = new Vector3();
    startingIndex = new int(0);
    selectedModel = (int)(m_scene->numModels() * Random::common().uniform());

	if (bGenerateVolumeSamples)
	{
		G3D::Vector3 sceneExtent = m_scene->m_maxBound - m_scene->m_minBound;
		
		G3D::Vector3 pt = G3D::Vector3(m_random.uniform(m_scene->m_minBound.x, m_scene->m_maxBound.x),
									   m_random.uniform(m_scene->m_minBound.y, m_scene->m_maxBound.y),
									   m_random.uniform(m_scene->m_minBound.z, m_scene->m_maxBound.z));

		while (bPreventOOBDisplacement && pointInsideEntity(pt))
		{
			pt = G3D::Vector3(m_random.uniform(m_scene->m_minBound.x, m_scene->m_maxBound.x),
							  m_random.uniform(m_scene->m_minBound.y, m_scene->m_maxBound.y),
							  m_random.uniform(m_scene->m_minBound.z, m_scene->m_maxBound.z));
		}

		*P = pt;
	}
	else
	{
		Array<Vector3> vertices = getRandomPoint(selectedModel, P, N, baryWeights, startingIndex);
	}
	SceneSample ss = SceneSample(*P, *N);

	delete baryWeights;
	delete P;
	delete N;
	delete startingIndex;

	return ss;
}

void App::clearAllActors(){
	if (actors.size() > 0)
	{
		//removeWidget(actors[0].getManipulator());
		actors.clear();
		randomPoints.clear();
		//resetParticles();
		extrapolationT = 0;
	}
}

void App::startOptimizationPasses()
{
	if (!sampleSetLoaded() && !probeStructureLoaded())
	{
		return;
	}

	sampleSet->outputToLog = logSampleSet;
	numPassesLeft = std::atoi(tbNumPassesLeft.c_str());	

	std::string settingsFilePath = "../data-files/scripts/optimizationSettings.txt";

	currentOptimization.bWaitingForRenderingFinished = false;
	currentOptimization.lastRenderEndTime = getFileLastModifiedTime(settingsFilePath.c_str());
}

float App::computeError(bool outputToLog)
{
	std::fstream currentFile, referenceFile, logFile;
    String errorLogPath = currentOptimizationFolderPath() + "/errorlog.txt";

	currentFile.open(   (currentOptimizationFolderPath() + "/values.txt"    ).c_str(), std::fstream::in);
	referenceFile.open( (currentOptimizationFolderPath() + "/ref_values.txt").c_str(), std::fstream::in);
    logFile.open(       errorLogPath.c_str(),                                          std::fstream::out | std::fstream::app);

	float error = 0;

    std::string currentLine, referenceLine;
    while (std::getline(currentFile, currentLine))
	{
		std::getline(referenceFile, referenceLine);
		float currentValue = std::stof(currentLine.c_str());
		float referenceValue = std::stof(referenceLine.c_str());
	
		error += pow(currentValue - referenceValue, 2);
	}

	currentFile.close();
	referenceFile.close();

    if (!logFile.is_open())
    {
        debugPrintf("Couldn't open error log file at %s\n", errorLogPath);
    }

	if (outputToLog)
	{
		logFile.precision(20);
		logFile << error << std::endl;
	}

	logFile.close();


	return error;

}


void App::createTempProbeStructure(G3D::Array<G3D::Vector3>& probePositions)
{
	ProbeStructure tempStructure = ProbeStructure(m_scene->m_name, "temp", probePositions.size(), EProbeStructureType::WeightedNearestNeighbour);

	for (int i = 0; i < probePositions.size(); ++i)
	{
		Probe* p = tempStructure.getProbe(i);

		p->frame.translation = probePositions[i];
	}
	tempStructure.updateProbes(true);

}

bool App::displacementCrossesSurface(G3D::Vector3 startPoint, G3D::Vector3 displacement, TriTree::Hit& hit)
{
	Ray ray = Ray::fromOriginAndDirection(startPoint, displacement.unit(), 0.0f, 100.f);

	int score = 0;
	int counter = 0;
	if (m_triTree.intersectRay(ray, hit, TriTree::DO_NOT_CULL_BACKFACES))
	{
		return (hit.distance <= displacement.length());
	}
	return false;
}

bool App::displacementCrossesSurface(G3D::Vector3 startPoint, G3D::Vector3 displacement)
{
	TriTree::Hit hit;
	return displacementCrossesSurface(startPoint, displacement, hit);

}

bool App::pointInsideEntity(G3D::Vector3 point, TriTree::Hit& hit)
{
	Ray ray = Ray::fromOriginAndDirection(point, Vector3(1,0,0), 0.0f, 100.f);

	int score = 0;
	int counter = 0;
	TriTree::Hit h;
	while (m_triTree.intersectRay(ray, h, TriTree::DO_NOT_CULL_BACKFACES))
	{
		if (counter++ == 0)
		{
			hit = h;
		}
		if (hit.backface)
		{
			score -= 1;
		}
		else
		{
			score += 1;
		}
		

		ray = Ray::fromOriginAndDirection(ray.origin() + (hit.distance + 0.01) * ray.direction(), ray.direction());
	}

	return score <= 0;
}

bool App::pointInsideEntity(G3D::Vector3 point)
{
	TriTree::Hit hit;
	return pointInsideEntity(point, hit);
}

G3D::Array<G3D::Vector3> App::generateRandomPositions(int NumberOfPositions)
{
	G3D::Array<G3D::Vector3> toReturn;
	float maxDistance = 0.5;
	float offset = 0.5;
	G3D::Vector3 min = m_scene->m_minBound + Vector3(offset, offset, offset);
	G3D::Vector3 max = m_scene->m_maxBound - Vector3(offset, offset, offset);

	for (int i = 0; i < NumberOfPositions; ++i)
	{
		while(true)
		{
			bool ok = true;
			G3D::Vector3 possiblePos = G3D::Vector3(Random::common().uniform(min.x, max.x),
									   				Random::common().uniform(min.y, max.y),
									   				Random::common().uniform(min.z, max.z));

			for (int j = 0; j < toReturn.size(); ++j)
			{
				float distance = (possiblePos - toReturn[j]).length();

				if (distance < maxDistance)
				{
					ok = false;
				}
			}


			if (pointInsideEntity(possiblePos))
			{
				ok = false;
			}

			if (ok)
			{
				toReturn.append(possiblePos);
				break;
			}
		}
	}

	return toReturn;
}


void App::findBestInitialConditions()
{
	debugPrintf("Attempting to find minimalest error...\n");

	int NumberOfProbes = std::atoi(m_sNumICProbes.c_str());
	int NumberOfTries = std::atoi(m_sNumICTries.c_str());

	float bestError = 99999;
	G3D::Array<G3D::Vector3> bestPositions;
	bestPositions.resize(NumberOfProbes);

	for (int tryNumber = 0; tryNumber < NumberOfTries; ++tryNumber)
	{
		G3D::Array<G3D::Vector3> positions = generateRandomPositions(NumberOfProbes);


		for (G3D::Vector3& v : positions)
		{
			m_probeStructure->addProbe(v);
		}
		m_probeStructure->updateAll(bShowOptimizationOutput);

		computeSamplesRGB();

		float error = computeError(false);
		debugPrintf("- Attempt #%d, error : %f", tryNumber, error);

		if (error < bestError)
		{
			debugPrintf(" (best yet!)");
			bestError = error;
			bestPositions = positions;
		}
		debugPrintf("\n");

		for (int i = 0; i < NumberOfProbes; ++i)
		{
			m_probeStructure->removeProbe(m_probeStructure->probeCount() - 1);
		}

	}
	debugPrintf("error: %f \n", bestError);

	for (G3D::Vector3& v : bestPositions)
	{
		debugPrintf("v: %s \n", v.toString().c_str());
		m_probeStructure->addProbe(v);

	}
	m_probeStructure->updateAll(bShowOptimizationOutput);

}

void App::logOptimizationIteration(float error)
{
	String path = currentOptimizationFolderPath() + "/infos.txt";

	std::fstream infoFile;
	infoFile.open(path.c_str(), std::fstream::out | std::fstream::app);

	if (infoFile.is_open())
	{
		if (currentOptimization.iteration == 0)
		{		
			infoFile << "optimizing : " << m_probeStructure->m_name.c_str() << std::endl << std::endl;
		}

		infoFile << "iteration " << currentOptimization.iteration << std::endl;
		infoFile << "error " << error << std::endl;

		for (int i = 0; i < m_probeStructure->probeList.size(); ++i)
		{
			infoFile << "Probe " << i << " : " << m_probeStructure->probeList[i]->getPosition() << std::endl;
		}
	}
	infoFile << std::endl;
	infoFile.close();
}

std::vector<float> App::tryOptimization()
{
	G3D::StopWatch sw("OuterOptimization");
	sw.setEnabled(true);

	sw.after("Starting optimization pass #" + String(currentOptimization.iteration));
	G3D::String sceneName = scenePane.selectedSceneList->selectedValue();
	G3D::String probeStructureName = scenePane.probeStructureList->selectedValue();

	int numSamples = std::atoi(numOptimizationSamples.c_str());
	int numCoeffs = std::atoi(optimizationSHBand.c_str());
	std::stringstream ss;
	std::vector<float> displacements;
	sampleSet->oneRowPerSHBand = bOneRowPerSHBand;

	if (bOptimizeForCoeffs)
	{
		sampleSet->generateInterpolatedCoefficientsFromProbes(numSamples, numCoeffs, currentOptimizationFolderPath() + "/values.txt", 0);
	}
	else
	{
		sampleSet->generateRGBValuesFromProbes(numSamples, numCoeffs, currentOptimizationFolderPath() + "/values.txt", 0);
	}

	sw.after("Generated RGB values from probes");

    float error = computeError(false);
	logOptimizationIteration(error);
	sw.after("Error : " + String(error));

	currentOptimization.iteration++;

	if (currentOptimization.errors.size())
	{
		if (error > currentOptimization.errors[currentOptimization.errors.size() - 1])
		{
			currentOptimization.consecutiveFailures++;
		}
		else
		{
			currentOptimization.consecutiveFailures = 0;
		}
	}
	
	if (bPreventErrorIncrease && currentOptimization.consecutiveFailures >= 1)
	{
		return displacements;
	}
	computeError(true);

	currentOptimization.errors.push_back(error);

	sampleSet->coeffReference = bOptimizeForCoeffs;
	displacements = sampleSet->tryOptimizationPass(numSamples, numCoeffs, bOptimizeWithMitsubaSamples, currentOptimizationFolderPath());
	sw.after("Finished optimization pass");
	return displacements;

}

void App::computeSamplesRGB()
{
	int numSamples = std::atoi(numOptimizationSamples.c_str());
    int numCoeffs = std::atoi(optimizationSHBand.c_str());
    String outputFile = currentOptimizationFolderPath() + "/values.txt";
	sampleSet->oneRowPerSHBand = bOneRowPerSHBand;

	if (bOptimizeForCoeffs)
	{
		sampleSet->generateInterpolatedCoefficientsFromProbes(numSamples, numCoeffs, outputFile, 0);
	}
	else
	{
		sampleSet->generateRGBValuesFromProbes(numSamples, numCoeffs, outputFile, 0);
	}
}
void App::computeSamplesRGBRef()
{
	if (!sampleSetLoaded())
	{
		debugPrintf("No sample set loaded");
		return;
	}

	int numSamples = std::atoi(numOptimizationSamples.c_str());
	int numCoeffs = std::atoi(optimizationSHBand.c_str());
	String outputFile = currentOptimizationFolderPath() + "/ref_values.txt";
	sampleSet->oneRowPerSHBand = bOneRowPerSHBand;
	if (bOptimizeWithMitsubaSamples)
	{
		if (bOptimizeForCoeffs)
		{
			copyFile(currentSampleSetPath() + "/InterpolatedCoeffs.txt", outputFile);
		}
		else
		{
			sampleSet->generateRGBValuesFromSamples(numSamples, outputFile, 0);
		}
	}
	else
	{
		if (bOptimizeForCoeffs)
		{
			sampleSet->generateInterpolatedCoefficientsFromProbes(numSamples, numCoeffs, outputFile, 0);

		}
		else
		{
			sampleSet->generateRGBValuesFromProbes(numSamples, numCoeffs, outputFile, 0);
		}
	}

	std::fstream infoFile((currentOptimizationFolderPath() + "/infos.txt").c_str(), std::fstream::out | std::fstream::app);
	infoFile << "refProbeStructure " << m_probeStructure->name().c_str() << std::endl;
	infoFile << "numSamples " << numSamples << std::endl << std::endl;
	infoFile.close();

	bTakeRefScreenshot = true;
}

void App::computeTriplets()
{
	int numSamples = std::atoi(numOptimizationSamples.c_str());
    int numCoeffs = std::atoi(optimizationSHBand.c_str());
    sampleSet->generateTriplets(numSamples, numCoeffs, currentOptimizationFolderPath() + "triplets.txt", 0, false);
}

void App::addOneActor()
{
	addActor("bunny", sphereModel /*sceneModel*/, 
			 Vector3( std::stof(actorSpawnX.c_str()), std::stof(actorSpawnY.c_str()), std::stof(actorSpawnZ.c_str()) ),
			 0.1f, shared_ptr<Texture>(), false, Vector3(1,1,1));
}

void App::addOneActorSq()
{
	for (int i = 0; i < m_scene->m_models.size(); ++i)
	{
		addActor("square", m_scene->m_models[i], 
				 Vector3(std::stof(actorSpawnX.c_str()), std::stof(actorSpawnY.c_str()), std::stof(actorSpawnZ.c_str())),
				 1.0f, shared_ptr<Texture>(), false, m_scene->m_colors[i]);
	}
	//addActor("square", sceneModel, 
	//		 Vector3(stof(actorSpawnX.c_str()),stof(actorSpawnY.c_str()),stof(actorSpawnZ.c_str())), 
	//		 1.0f, shared_ptr<Texture>());
}

void App::displaceProbes()
{
	//G3D::Vector3 displacementV = Vector3(0.1, 0, 0).unit();

	//std::vector<float> displacement;
	//displacement.push_back(displacementV.x);
	//displacement.push_back(displacementV.y);
	//displacement.push_back(displacementV.z);

	//m_probeStructure->displaceProbesWithGradient(displacement, 0.1f);
	//m_probeStructure->saveCoefficients();

	//return;
	std::vector<G3D::Vector3> originalPositions;

	for (int i = 0 ; i < m_probeStructure->probeCount(); ++i)
	{
		originalPositions.push_back(m_probeStructure->getProbe(i)->frame.translation);
	}


	int numTries = 100;
	for (int i = 1; i <= numTries; ++i)
	{

		for (int j = 0 ; j < m_probeStructure->probeCount(); ++j)
		{
			m_probeStructure->getProbe(j)->setPosition(originalPositions[j]);
		}
		//m_probeStructure->updateProbes(true);
        m_probeStructure->savePositions(false);
        m_probeStructure->generateProbes("all", false, true, bShowOptimizationOutput);
        m_probeStructure->extractSHCoeffs(true, true);
		//updateProbeStructure();

		G3D::Vector3 displacementV = Vector3(1,1,1).unit();

		std::vector<float> displacement;
		displacement.push_back((float)i * 0.001f * displacementV.x);
		displacement.push_back((float)i * 0.001f * displacementV.y);
		displacement.push_back((float)i * 0.001f * displacementV.z);

		Array<G3D::Vector3> renderedCoeffs, optimCoeffs;

		m_probeStructure->displaceProbesWithGradient(displacement, 0.1f);

        m_probeStructure->savePositions(false);
		//computeSamplesRGB();
		//computeError("C:/temp/errorlog.txt");
		optimCoeffs = m_probeStructure->getProbe(0)->coeffs;

		//m_probeStructure->updateProbes(true);
		//updateProbeStructure();
        m_probeStructure->generateProbes("all", false, true, bShowOptimizationOutput);
        m_probeStructure->extractSHCoeffs(true, true);
		//computeSamplesRGB();
		//computeError("C:/temp/errorlog2.txt");
		renderedCoeffs = m_probeStructure->getProbe(0)->coeffs;

		float err = 0;
		for (int j = 0; j < optimCoeffs.size(); ++j)
		{
			G3D::Vector3 diff = optimCoeffs[j] - renderedCoeffs[j];

			err += diff.dot(diff);


		}
		std::fstream outFile;
		outFile.open("C:/temp/errorlog.txt", std::fstream::out | std::fstream::app);
		outFile << err << "\n";
		outFile.close();
	}
}


