#include "App.h"
#include "SceneSampleSet.h"
#define DEAV3(x) debugPrintf(#x); for(int num = 0; num < x.size(); num++){ debugPrintf(", [num]: (%s)\n",x[num].toString().c_str()); debugPrintf("\n");}



Vector3 App::getRandomPointInScene(){
	float x = r.uniform(-5,5);
	float y = r.uniform(0,5);
	float z = r.uniform(-5,5);

	return Vector3(x,y,z);
}

Array<Vector3> App::getRandomPoint(int modelNumber, Vector3* P, Vector3* N, Vector3* barycentricWeights, int* s){
    Array<Vector3> toReturn = Array<Vector3>();
	//Random r(System::time());

	MeshShape* meshShape = m_scene.getMeshShape(modelNumber);

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

void App::sampleScenePoint(int* _selectedModel, Vector3 *_P, Vector2* _UV, Vector3* _N, int sampleID){

    Vector3* baryWeights = new Vector3();
    Vector3* P = new Vector3();
    Vector3* N = new Vector3();
    startingIndex = new int(0);
    selectedModel = (int)(m_scene.numModels() * r.uniform());
    //selectedModel = 3;
    Array<Vector3> vertices = getRandomPoint(selectedModel, P, N, baryWeights, startingIndex);


    *_selectedModel = selectedModel;
    //*_UV = UV;
    *_N = *N;
	*_P = *P;

		


	Vector2 IJ = Vector2(0, 0);
    Color4 uvColor = Color4(0, 0, 0, 0);


	SceneSample ss = SceneSample(*P, *N);
	sampleSet->addSample(ss);
	//if (saveSample) ss->writeToFile("C:/libraries/g3d/samples/aTest/data-files/Scenes/" + 
	//								scenePane.selectedSceneList->selectedValue() + "/SampleSets/" + 
	//								*SampleSetOutputName + "/SamplePositions.txt");

	delete baryWeights;
	delete P;
	delete N;
	delete startingIndex;
}

void App::clearAllActors(){
	if (actors.size() > 0)
	{
		removeWidget(actors[0].getManipulator());
		actors.fastClear();
		randomPoints.clear();
		//resetParticles();
		extrapolationT = 0;
	}
}

void App::startOptimizationPasses()
{
	if (sampleSet)
	{
		std::fstream infoFile;
		infoFile.open("C:/temp/CurrentOptimization/InitialConditions.txt", std::fstream::out | std::fstream::app | std::fstream::in);
		
		if (!infoFile)
		{
			infoFile.open("C:/temp/CurrentOptimization/InitialConditions.txt", std::fstream::in | std::fstream::out | std::fstream::trunc);
			infoFile << m_scene.m_name.c_str() << std::endl;
			infoFile << probeStructure->m_name.c_str() << std::endl;

			for (int i = 0; i <  probeStructure->probeList.size(); ++i)
			{
				infoFile << "Probe " << i << " : " << probeStructure->probeList[i]->getPosition() << std::endl;
			}
		}
		sampleSet->outputToLog = logSampleSet;
		numPassesLeft = std::atoi(tbNumPassesLeft.c_str());
	}
	
}
float App::computeError(std::string logFilePath)
{
	std::fstream currentFile, referenceFile;

	currentFile.open("C:/temp/samplesRGB.txt", std::fstream::in);
	referenceFile.open("C:/temp/samplesRGB_ref.txt", std::fstream::in);

	std::string currentLine, referenceLine;

	float error = 0;

	while (std::getline(currentFile, currentLine))
	{
		std::getline(referenceFile, referenceLine);
		float currentValue = std::stof(currentLine.c_str());
		float referenceValue = std::stof(referenceLine.c_str());
	
		error += pow(currentValue - referenceValue, 2);
	}

	currentFile.close();
	referenceFile.close();

	std::fstream logFile;
	logFile.open(logFilePath.c_str(), std::fstream::out | std::fstream::app);
	if (!logFile)
	{
		// create file because it does not exist
		logFile.open(logFilePath.c_str(),  std::fstream::in | std::fstream::out | std::fstream::trunc);
	}
	logFile.precision(20);
	logFile << error << std::endl;
	logFile.close();


	return error;

}


void App::createTempProbeStructure(G3D::Array<G3D::Vector3>& probePositions)
{
	ProbeStructure tempStructure = ProbeStructure(m_scene.m_name, "temp", probePositions.size(), EProbeStructureType::WeightedNearestNeighbour);

	for (int i = 0; i < probePositions.size(); ++i)
	{
		Probe* p = tempStructure.getProbe(i);

		p->frame.translation = probePositions[i];
	}
	tempStructure.updateProbes(true);

}

G3D::Array<G3D::Vector3> App::generateRandomPositions(int NumberOfPositions)
{
	G3D::Array<G3D::Vector3> toReturn;
	float maxDistance = 4;
	float offset = 0.5;
	G3D::Vector3& min = m_scene.m_minBound + Vector3(offset, offset, offset);
	G3D::Vector3& max = m_scene.m_maxBound - Vector3(offset, offset, offset);

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
	int NumberOfProbes = std::atoi(m_sNumICProbes.c_str());
	int NumberOfTries = std::atoi(m_sNumICTries.c_str());

	float bestError = 99999;
	G3D::Array<G3D::Vector3> bestPositions;
	bestPositions.resize(NumberOfProbes);

	for (int tryNumber = 0; tryNumber < NumberOfTries; ++tryNumber)
	{
		G3D::Array<G3D::Vector3> positions = generateRandomPositions(NumberOfProbes);

		createTempProbeStructure(positions);

		initializeProbeStructure(m_scene.m_name, "temp");
		
		computeSamplesRGB();

		float error = computeError("C:/temp/errorlog.txt");

		if (error < bestError)
		{
			bestError = error;
			bestPositions = positions;
		}

		std::fstream positionsFile;
		positionsFile.open("C:/temp/positionsLog/" + std::to_string(tryNumber) + ".txt", std::fstream::out);
		for (int i = 0; i < NumberOfProbes; ++i)
		{
			positionsFile << positions[i] << std::endl;
		}
		positionsFile.close();
	}
	debugPrintf("error: %f \n", bestError);

	for (G3D::Vector3& v : bestPositions)
	{
		debugPrintf("v: %s \n", v.toString().c_str());

	}
	// generate 

	createTempProbeStructure(bestPositions);
	//exit(0);
	initializeProbeStructure(m_scene.m_name, "temp");
}


void App::tryOptimization()
{
	G3D::StopWatch sw;
	G3D::String sceneName = scenePane.selectedSceneList->selectedValue();
	G3D::String probeStructureName = scenePane.probeStructureList->selectedValue();

	int numSamples = std::atoi((*samplesToSave).c_str());
	//sampleSet->generateRGBValuesFromProbes(numSamples, false, 0);

	std::stringstream ss;
	std::vector<float> displacements;
	if (useMatlabOptimization)
	{
		computeSamplesRGB();
		computeTriplets();
		ss << "set PATH=%PATH%;C:\\Program Files\\MATLAB\\MATLAB Runtime\\v85\\runtime\\win64";
		ss << " && C:\\libraries\\sparsematrixfitting_onecolumn\\for_testing\\sparsematrixfitting_onecolumn.exe";
		system(ss.str().c_str());

		std::fstream file;
		file.open("C:/temp/dp.txt", std::fstream::in);
		std::string line;
		while(std::getline(file, line))
		{
			displacements.push_back(std::stof(line.c_str()));
		}
	}
	else
	{
		displacements = sampleSet->tryOptimizationPass(numSamples, false);
	}
	sw.after("Performed optimization step");
	if (displacements.size() > 0)
	{
		probeStructure->displaceProbesWithGradient(displacements);
		//sampleSet->generateRGBValuesFromProbes(numSamples, false,0);
		//float error = computeError("C:/temp/CurrentOptimization/errorlog.txt");
		//sw.after("Computed error");

		//probeStructure->applyOffsetToProbes(displacements);
		//sw.after("Applied offset to probes");
		//probeStructure->updateProbes(true);
		//sw.after("Computed new probe textures and coefficients");


		//updateProbeStructure();
		//sw.after("Reloaded probe structure");
		//updateSampleSet();
	}
}

void App::computeSamplesRGB()
{
	int numSamples = std::atoi((*samplesToSave).c_str());
	sampleSet->generateRGBValuesFromProbes(numSamples, false, 0);
}
void App::computeSamplesRGBRef()
{
	int numSamples = std::atoi((*samplesToSave).c_str());
	sampleSet->generateRGBValuesFromProbes(numSamples, true, 0);
}

void App::computeTriplets()
{
	int numSamples = std::atoi((*samplesToSave).c_str());
	sampleSet->generateTriplets(numSamples, 0);
}


void App::computeSceneSamples(){
    // generate n random points and store their colors, compare to the color in the corresponding texture , add to a cost function
	int n = atoi(samplesToSave->c_str());
    Vector2 UV, IJ;
    Vector3 N, P;
    int selectedModel;

    for (int i = 0; i < n; i++)
	{
		clearAllActors();
        sampleScenePoint(&selectedModel, &P, &UV, &N, i);
    }
	sampleSet->save();
}




void App::addOneActor()
{
	addActor("bunny", bunnyModel /*sceneModel*/, 
			 Vector3( std::stof(actorSpawnX.c_str()), std::stof(actorSpawnY.c_str()), std::stof(actorSpawnZ.c_str()) ),
			 0.1f, shared_ptr<Texture>(), true, Vector3(1,1,1));
}

void App::addOneActorSq()
{
	for (int i = 0; i < m_scene.m_models.size(); ++i)
	{
		addActor("square", m_scene.m_models[i], 
				 Vector3(std::stof(actorSpawnX.c_str()), std::stof(actorSpawnY.c_str()), std::stof(actorSpawnZ.c_str())),
				 1.0f, shared_ptr<Texture>(), false, m_scene.m_colors[i]);
	}
	//addActor("square", sceneModel, 
	//		 Vector3(stof(actorSpawnX.c_str()),stof(actorSpawnY.c_str()),stof(actorSpawnZ.c_str())), 
	//		 1.0f, shared_ptr<Texture>());
}

void App::displaceProbes()
{

	std::vector<G3D::Vector3> originalPositions;

	for (int i = 0 ; i < probeStructure->probeCount(); ++i)
	{
		originalPositions.push_back(probeStructure->getProbe(i)->frame.translation);
	}


	int numTries = 100;
	for (int i = 1; i <= numTries; ++i)
	{

		for (int i = 0 ; i < probeStructure->probeCount(); ++i)
		{
			probeStructure->getProbe(i)->frame.translation = originalPositions[i];
		}
		probeStructure->updateProbes(true);
		updateProbeStructure();

		G3D::Vector3 displacementV = Vector3(1,1,1).unit();

		std::vector<float> displacement;
		displacement.push_back((float)i * 0.001f * displacementV.x);
		displacement.push_back((float)i * 0.001f * displacementV.y);
		displacement.push_back((float)i * 0.001f * displacementV.z);

		Array<G3D::Vector3> renderedCoeffs, optimCoeffs;

		probeStructure->displaceProbesWithGradient(displacement);
		computeSamplesRGB();
		//computeError("C:/temp/errorlog.txt");
		optimCoeffs = probeStructure->getProbe(0)->coeffs;

		probeStructure->updateProbes(true);
		updateProbeStructure();
		computeSamplesRGB();
		//computeError("C:/temp/errorlog2.txt");
		renderedCoeffs = probeStructure->getProbe(0)->coeffs;

		float err = 0;
		for (int i = 0; i < optimCoeffs.size(); ++i)
		{
			G3D::Vector3 diff = optimCoeffs[i] - renderedCoeffs[i];

			err += diff.dot(diff);


		}
		std::fstream outFile;
		outFile.open("C:/temp/errorlog.txt", std::fstream::out | std::fstream::app);
		outFile << err << "\n";
		outFile.close();
	}
}


