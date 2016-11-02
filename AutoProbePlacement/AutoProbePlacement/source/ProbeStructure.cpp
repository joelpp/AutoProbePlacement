#include <sstream>

#include "ProbeStructure.h"
#include "App.h"
#include "Helpers.h"
#include "SH.h"

App* App::instance;
G3D::Array<G3D::String> ProbeStructure::typeMap;

inline void ProbeStructure::createTypeMap()
{
	typeMap.resize(EProbeStructureType::NUM_TYPES);

	typeMap[EProbeStructureType::Trilinear] = "trilinear";;
	typeMap[EProbeStructureType::Closest] = "closest";
	typeMap[EProbeStructureType::Tetrahedral] = "tetrahedral";
	typeMap[EProbeStructureType::WeightedNearestNeighbour] = "wNN";
}

ProbeStructure::ProbeStructure()
{

}

ProbeStructure::ProbeStructure(String sceneName, String probeStructureName, int numProbes, EProbeStructureType type)
{
	createTypeMap();

	this->m_sceneName = sceneName;
	this->m_name = probeStructureName;
	this->probeStructurePath = "../data-files/Scenes/" + sceneName + "/ProbeStructures/" + probeStructureName;;
	this->m_dimensions.push_back(numProbes);


	for (int i = 0; i < numProbes; ++i)
	{
		this->probeList.push_back(new Probe(i, probeStructurePath));
		this->probeList[i]->bNeedsUpdate = true;
	}
	this->m_type = type;

	std::fstream infoFile;
	infoFile.open((probeStructurePath + "/info.txt").c_str(), std::fstream::out);
	infoFile << "type wNN" << std::endl;
	infoFile << "dimensions " << numProbes << std::endl;
	infoFile.close();


}

ProbeStructure::ProbeStructure(String sceneName, String probeStructureName)
{
	createTypeMap();

	this->m_sceneName = sceneName;
	this->m_name = probeStructureName;
	this->probeStructurePath = "../data-files/Scenes/" + sceneName + "/ProbeStructures/" + probeStructureName;;

	loadProbeStructureInfo();
	makeProbeList();

    m_NumCoefficients = 9;
    m_NumColors = 3;

	if (m_type == EProbeStructureType::Tetrahedral)
	{
		makeLineArray();
		makeTetrahedraArray();
		loadTetrahedraNeighbors();
		loadHullFaces();
		computeTetFaceNormals();
		computeTetVertexNormals();
	}
 
    uploadToGPU();
}


void ProbeStructure::loadProbeStructureInfo()
{
	String probeInfoPath = probeStructurePath + "/info.txt";
	std::fstream structureInformationFile;

    structureInformationFile.open(probeInfoPath.c_str(), std::ios::in);
	std::string line;

	while ( std::getline(structureInformationFile, line) )
	{
		Array<String> splitLine = stringSplit(String(line.c_str()), ' ');
		String param = splitLine[0]; 
		splitLine.remove(0);
		if (param == "type")
		{
			
			String& s = splitLine[0];

			this->m_type = EProbeStructureType((int)typeMap.findIndex(s) );

		}
		if (param == "step")
		{
			this->m_step = std::stof(splitLine[0].c_str());
		}
		else if (param == "dimensions")
		{
			for (int i = 0; i < splitLine.size(); ++i)
			{
				this->m_dimensions.push_back(std::stoi(splitLine[i].c_str()));
			}
		}
		else if (param == "firstProbePosition")
		{
			this->m_firstProbePosition[0] = std::stof(splitLine[0].c_str());
			this->m_firstProbePosition[1] = std::stof(splitLine[1].c_str());
			this->m_firstProbePosition[2] = std::stof(splitLine[2].c_str());
		}
        else if (param == "gamma")
        {
            this->m_gamma = std::stof(splitLine[0].c_str());
        }
        else if (param == "sampleCount")
       { 
            this->m_NumSamples = std::stoi(splitLine[0].c_str());
        }
        else if (param == "width")
        {
            this->m_width = std::stoi(splitLine[0].c_str());
        }
        else if (param == "height")
        {
            this->m_height = std::stoi(splitLine[0].c_str());
        }
        else if (param == "integrator")
		{
			this->m_integrator = splitLine[0];
		}



	}
}

std::fstream openFile(String& path)
{
	std::fstream returnFile;

	returnFile.open(path.c_str(), std::ios::in);
	if (!returnFile.is_open())
	{
		debugPrintf("Could not open file: %s", path.c_str());
		exit(1);
	}
	return returnFile;
}
void ProbeStructure::makeProbeList()
{
    //debugPrintf("Making probe list...\n");
	probeList = Array<Probe*>();
    bool loadCoefficients = true;

	std::fstream probeInfo = openFile(probeStructurePath + "/probeList.txt");

    std::string line, gradientLine;
    Array<Point3> positions = Array<Point3>();
    int probeCount = 0;

    while ( std::getline(probeInfo,line) )
	{
        Array<String> coords = stringSplit(String(line.c_str()), ' ');

		if (coords[0] == "#")
		{
			continue;
		}
        positions.append(Point3(std::stof(coords[0].c_str()), std::stof(coords[1].c_str()), std::stof(coords[2].c_str())));

        probeCount++;
    }

    probeInfo.close();
    Probe* aTestProbe; 

    for (int i = 0; i < probeCount; i++)
	{

        aTestProbe = new Probe(i, probeStructurePath);
		
        String filename;
                
        char str[16];
        sprintf(str, "%d", i);
        aTestProbe->position = positions[i];
                
        if (loadCoefficients)
		{
			////////////////////////////////////
			//////// LOAD COEFFICIENTS /////////
			////////////////////////////////////
			filename = probeStructurePath + "/coefficients/Coefficients_";
            filename.append(str);
            filename.append(".txt");

            std::string line;
            std::string line2;
            std::string line3;
            std::fstream coeffFile;
            coeffFile.open(filename.c_str(), std::ios::in);

			aTestProbe->coeffs = Array<Vector3>();

            while ( std::getline(coeffFile,line) )
			{
                std::getline(coeffFile, line2);
                std::getline(coeffFile, line3);
                Vector3 toAdd = Vector3(std::stof(line.c_str()), std::stof(line2.c_str()), std::stof(line3.c_str()));
                aTestProbe->coeffs.append(toAdd);
            }
            coeffFile.close();

			////////////////////////////////////
			//////// LOAD GRADIENTS ////////////
			////////////////////////////////////

			filename = probeStructurePath + "/coefficients/CoefficientsGradients_";
			filename.append(str);
			filename.append(".txt");

			aTestProbe->coeffGradients = Array< Array<Vector3> >();
			Array<Vector3> temporaryArray;
			coeffFile.open(filename.c_str(), std::ios::in);

			if (!coeffFile.is_open())
			{
				for (int j = 0; j < 3; ++j)
				{
					temporaryArray.push_back(Vector3(0,0,0));
				}
				for (int i = 0; i < 9; ++i)
				{
					aTestProbe->coeffGradients.push_back(temporaryArray);
				}
			}
			else
			{
				int counter = 0;
				while (std::getline(coeffFile, line))
				{
					while (line[0] == '/')
					{
						std::getline(coeffFile, line);
					}

					std::getline(coeffFile, line2);
					std::getline(coeffFile, line3);
					temporaryArray.push_back(Vector3(std::stof(line.c_str()), std::stof(line2.c_str()), std::stof(line3.c_str())));
					//temporaryArray.push_back(Vector3(Random::common().uniform(-0.15f, 0.15f), Random::common().uniform(-0.15f, 0.15f), Random::common().uniform(-0.15f, 0.15f)));

					if (++counter == 3)
					{
						aTestProbe->coeffGradients.push_back(temporaryArray);
						temporaryArray.clear();
						counter = 0;
					}

				}
				coeffFile.close();
			}
        }

		aTestProbe->m_sphere = Sphere(positions[i], 0.1f);

		
		aTestProbe->frame = CFrame::fromXYZYPRDegrees(aTestProbe->position.x, aTestProbe->position.y, aTestProbe->position.z, 0, 0, 0);


		aTestProbe->manipulator = ProbeManipulator::create();
		aTestProbe->manipulator->setFrame(aTestProbe->frame);
		aTestProbe->manipulator->setControlFrame(aTestProbe->frame);
		aTestProbe->manipulator->p = aTestProbe;

		aTestProbe->bNeedsUpdate = false;
        probeList.append(aTestProbe);
		probeMap.insert(std::pair<Vector3, Probe*>(aTestProbe->position, aTestProbe));
    }
}


void ProbeStructure::reset(){
	makeProbeList();
}

void ProbeStructure::makeLineArray(){
    lineArray = Array<LineSegment>();

	std::fstream lineIn = openFile(probeStructurePath + "/tetgenOutput/probeList.1.edge");

    std::string line; 
    std::getline(lineIn,line);

    Array<String> splitLine = stringSplit(String(line.c_str()), ' ');
    while (std::getline(lineIn, line)){
        bool readIndex = false;
        bool readEdge1 = false;
        bool readEdge2 = false;
        int edge1;
        int edge2;

        splitLine = stringSplit(String(line.c_str()), ' ');
        for (int i = 0; i < splitLine.size(); i++){

            int readNumber = atoi(splitLine[i].c_str());

            if (readNumber != 0){
                readNumber -= 1;

                if (!readIndex) readIndex = true;

                else if (!readEdge1){
                    readEdge1 = true;
                    edge1 = readNumber;
                    continue;
                }
                else if(!readEdge2){
                    readEdge2 = true;
                    edge2 = readNumber;
                    continue;
                }
                else continue;
            }
        }
        lineArray.append(LineSegment::fromTwoPoints(probeList[edge1]->position, probeList[edge2]->position));
    }
    lineIn.close();
    
}


void ProbeStructure::makeTetrahedraArray(){
    tetrahedraArray = Array<Tetrahedron*>();
	std::fstream lineIn = openFile(probeStructurePath + "/tetgenOutput/probeList.1.ele");

    std::string line; 
    Array<String> splitLine;	

    //first line useless
    std::getline(lineIn,line);


    int counter = 1;
    while (std::getline(lineIn, line)){
        bool readIndex = false;
        bool read1 = false;
        bool read2 = false;
        bool read3 = false;
        bool read4 = false;

        Array<int> myIndexes = Array<int>();

        splitLine = stringSplit(String(line.c_str()), ' ');
        if (splitLine[0].compare("#") == 0) break;

        //horrible
        for (int i = 0; i < splitLine.size(); i++){

            // -1 because tetgen starts numbering at 1
            int readNumber = atoi(splitLine[i].c_str());

            if (readNumber != 0){
                readNumber -= 1;

                if (!readIndex) readIndex = true;
                else if (!read1){
                    read1 = true;
                    myIndexes.append(readNumber);
                    continue;
                }
                else if(!read2){
                    read2 = true;
                    myIndexes.append(readNumber);
                    continue;
                }
                else if(!read3){
                    read3 = true;
                    myIndexes.append(readNumber);
                    continue;
                }
                else if(!read4){
                    read4 = true;
                    myIndexes.append(readNumber);
                    continue;
                }
                else continue;
            }

        }

        Tetrahedron* t = new Tetrahedron(counter, myIndexes);
        tetrahedraArray.append(t);
        counter++;

    }

    lineIn.close();
}


void ProbeStructure::loadTetrahedraNeighbors(){
    //debugPrintf("Loading tetrahedra neighbors\n");
    
    
	std::fstream lineIn = openFile(probeStructurePath + "/tetgenOutput/probeList.1.neigh");

    std::string line; 
    Array<String> splitLine;
    //first line useless

    std::getline(lineIn,line);
    int counter = 0;
    while (std::getline(lineIn, line)){
        bool readIndex = false;
        bool read1 = false;
        bool read2 = false;
        bool read3 = false;
        bool read4 = false;

        Array<int> myNeighbors = Array<int>();

        splitLine = stringSplit(String(line.c_str()), ' ');
        if (splitLine[0].compare("#") == 0) break;

        for (int i = 0; i < splitLine.size(); i++){

            // -1 because tetgen starts numbering at 1
            int readNumber = atoi(splitLine[i].c_str());
            // if (readNumber == -1) continue;
            if (readNumber != 0){
                if (readNumber != -1) readNumber -= 1;

                if (!readIndex) readIndex = true;
                else if (!read1){
                    read1 = true;
                    myNeighbors.append(readNumber);
                    continue;
                }
                else if(!read2){
                    read2 = true;
                    myNeighbors.append(readNumber);
                    continue;
                }
                else if(!read3){
                    read3 = true;
                    myNeighbors.append(readNumber);
                    continue;
                }
                else if(!read4){
                    read4 = true;
                    myNeighbors.append(readNumber);
                    continue;
                }
                else continue;
            }

        }

        tetrahedraArray[counter]->neighbors = myNeighbors;
        //debugPrintf("%s \n", tetrahedraArray[counter]->toString().c_str());
        counter++;
    }
    lineIn.close();

}

void ProbeStructure::loadHullFaces(){
    //debugPrintf("Loading Hull Faces\n");
    tetFaces = Array<Array<int> >();
	std::fstream lineIn = openFile(probeStructurePath + "/tetgenOutput/probeList.1.face");

    Array<String> splitLine;
    Array<int> goodSplit;
    std::string line;
    std::getline(lineIn,line);
    int counter = 0;
    while (std::getline(lineIn, line)){
        splitLine = stringSplit(String(line.c_str()), ' ');
        if (splitLine[0].compare("#") == 0) break;
        goodSplit = Array<int>();

        for (int i = 0; i < splitLine.size(); i++){
            if (splitLine[i].compare("") != 0) goodSplit.append(atoi(splitLine[i].c_str()));
        }

        if (goodSplit[4] != 1) continue;
        // debugPrintf(.c_str());
        // debugStringArray(splitLine);
        // debugIntArray(goodSplit);
        // debugPrintf("\n");
        tetFaces.append(goodSplit);
        counter++;
    }

    lineIn.close();
}


void ProbeStructure::computeTetFaceNormals(){
    //debugPrintf("Computing Tetrahedron face normals\n");

    tetFaceNormals = Array<Vector3>();
    tetFaceNormalPositions = Array<Vector3>();
    for (int i =0 ; i < tetFaces.size(); i++){
        Vector3 faceNormal;

        Array<int> vertices = Array<int>();
        vertices.append(tetFaces[i][1]-1);
        vertices.append(tetFaces[i][2]-1);
        vertices.append(tetFaces[i][3]-1);

        Vector3 line1 = probeList[vertices[0]]->position - probeList[vertices[1]]->position;
        Vector3 line2 = probeList[vertices[1]]->position - probeList[vertices[2]]->position;

        faceNormal = -line1.cross(line2).fastUnit();
        Vector3 midpoint = (1.0/3.0) * (probeList[vertices[0]]->position + probeList[vertices[1]]->position + probeList[vertices[2]]->position );

        tetFaceNormals.append(faceNormal);
        tetFaceNormalPositions.append(midpoint);
    }
}

void ProbeStructure::computeTetVertexNormals(){
    //debugPrintf("Computing tetrahedron vertex normals\n");

    //Iterate on probes
    for (int i = 0 ; i < probeList.size(); i++){
        //debugPrintf("Looking at probe #%d\n", i);
        Vector3 vertexNormal = Vector3(0,0,0);

        //Iterate on faces
        for (int j = 0; j < tetFaces.size(); j++){
            //debugPrintf("Looking at face %d\n",j);

            //Look at vertex indices forming the face
            for (int k = 1; k < 4; k++){

                //If our vertex index is on this face
                if (tetFaces[j][k] == i+1){
                    //debugPrintf("Match! This probe is tied to the face.\n");

                    vertexNormal += tetFaceNormals[j];
                }
            }
        }
        probeList[i]->normal = vertexNormal.fastUnit();
		//debugPrintf("probeList[%d].normal = %s\n", i, probeList[i].normal.toString().c_str());
    }
}

//std::string probeStructurePath = "RegularGrid";


void ProbeStructure::changeProbePosition(int probeIndex, Vector3 dp){
	//change probe position in list
	probeList[probeIndex]->position += dp;

	writeProbePositionsToFile();

	////re-render probe that has changed
	//char args[50];
	//sprintf(args, "%f %f %f %d", probeList[probeIndex].position.x, probeList[probeIndex].position.y, probeList[probeIndex].position.z, probeIndex);

	String command = String("cd C:/libraries/g3d/samples/aTest/data-files/scripts && python onecamera.py ");
	//command.append(args);
    //int result = system(command.c_str());
    int result = 0;

	//debugPrintf("result of onecamera: %d\n",result);

	////extract SH coeffs from new probe
	//char args2[10];
	//command = String("cd C:/libraries/g3d/samples/aTest/data-files/scripts && python extractprobecoefficients.py ");
	//sprintf(args2, "%d", probeIndex);
	//command.append(args2);
	//result = system(command.c_str());

	//debugPrintf("result of extractcoefficients: %d\n",result);


	//render new probe irr texture
	//command = String("cd C:/libraries/g3d/samples/aTest/data-files/scripts && python renderaprobe.py ");
	//command.append(args2);
	//result = system(command.c_str());
	//debugPrintf("result of renderaprobe: %d\n",result);

	//recalculate tetrahedrons
	writeProbePositionsToTetgenFile();
	command = String("cd C:/libraries/g3d/samples/aTest/data-files/tetgen/Debug && tetgen ../test/probeList.node -n -f -e");
	result = system(command.c_str());
	//debugPrintf("result of tetgen: %d\n",result);

	reset();

	//reload tetrahedrons
}

void ProbeStructure::changeProbeCoeff(int probeIndex, int coeffIndex, Vector3 amount){
	Probe *p = probeList[probeIndex];

	p->coeffs[coeffIndex] += amount;
}

bool ProbeStructure::writeProbePositionsToFile(){
	std::fstream posFile;
	posFile.open("../data-files/ProbeStructures/current/probeList.txt", std::fstream::out);
     for (int i = 0; i <probeList.size(); i++) posFile << probeList[i]->position.x << " " 
													   << probeList[i]->position.y << " " 
													   << probeList[i]->position.z << " " 
													   << "\n";
    posFile.close();
	return true;
}
bool ProbeStructure::writeProbePositionsToTetgenFile(){
	std::fstream posFile;
	posFile.open("../data-files/tetgen/current/probeList.node", std::fstream::out);
	posFile << probeList.size() << " 3 0 0\n";
    for (int i = 0; i <probeList.size(); i++) posFile << i+1 << " " 
													  << probeList[i]->position.x << " " 
													  << probeList[i]->position.y << " " 
													  << probeList[i]->position.z << " " 
													  << "\n";
    posFile.close();
	return true;
}

Probe* ProbeStructure::getProbe(int i)
{

	return probeList[i];

}

int ProbeStructure::probeCount(){
	return probeList.size();
}

Tetrahedron* ProbeStructure::getTetrahedron(int i){
	return tetrahedraArray[i];
}

Vector3 ProbeStructure::getFaceNormal(int i){
	return tetFaceNormals[i];
}

Vector3 ProbeStructure::getFaceNormalPosition(int i){
	return tetFaceNormalPositions[i];
}

LineSegment ProbeStructure::getLineSegment(int i){
	return lineArray[i];
}

int ProbeStructure::lineArraySize(){
	return lineArray.size();
}
int ProbeStructure::tetFaceNormalCount(){
	return tetFaceNormals.size();
}

void ProbeStructure::operator=(const ProbeStructure &p ){ 
	this->probeList = p.probeList;
}

int ProbeStructure::getProbeIndex(const G3D::Vector3& testPos)
{
	for (int i = 0; i < probeList.size(); ++i)
	{
		G3D::Vector3& pos = probeList[i]->position;
		if (pos == testPos)
		{
			return i;
		}
	}
	return -1;
}

G3D::Vector3 findProbe000(const G3D::Vector3& pos, float step)
{
	G3D::Vector3 toReturn;

	toReturn[0] = std::floorf((pos[0]) / step) * step;
	toReturn[1] = std::floorf((pos[1]) / step) * step;
	toReturn[2] = std::floorf((pos[2]) / step) * step;

	return toReturn;
}

G3D::Array<G3D::Vector3> ProbeStructure::getInterpolatingProbesCoords(const G3D::Vector3& pos, int step)
{
	G3D::Array<G3D::Vector3> toReturn;
	G3D::Vector3 probe000Pos = findProbe000(pos, (float) step);

	float fstep = (float)step;
	if (m_type == EProbeStructureType::Trilinear)
	{
		toReturn.push_back(probe000Pos);

		toReturn.push_back(probe000Pos + G3D::Vector3(fstep, 0, 0));
		toReturn.push_back(probe000Pos + G3D::Vector3(0, fstep, 0));
		toReturn.push_back(probe000Pos + G3D::Vector3(0, 0, fstep));
		toReturn.push_back(probe000Pos + G3D::Vector3(fstep, 0, fstep));
		toReturn.push_back(probe000Pos + G3D::Vector3(0, fstep, fstep));
		toReturn.push_back(probe000Pos + G3D::Vector3(fstep, fstep, 0));
		toReturn.push_back(probe000Pos + G3D::Vector3(fstep, fstep, fstep));
	}
	else if (m_type == EProbeStructureType::Trilinear)
	{
		G3D::Vector3 positionToReturn;
		float minDistance = 99999;

		for (int i = 0; i < probeCount(); ++i)
		{
			G3D::Vector3& probePosition = getProbe(i)->position;

			float distance = (pos - probePosition).length();

			if (distance < minDistance)
			{
				positionToReturn = probePosition;
				minDistance = distance;
				continue;
			}
		}

		toReturn.push_back(positionToReturn);
	}


	return toReturn;
}


G3D::Array<int> ProbeStructure::getInterpolatingProbeIndices(const G3D::Vector3& pos)
{
	G3D::Array<int> toReturn;
	if (m_type == EProbeStructureType::Closest)
	{
		int bestMatch = -1;
		float minDistance = 99999;

		for (int i = 0; i < probeCount(); ++i)
		{
			G3D::Vector3& probePosition = getProbe(i)->position;

			float distance = (pos - probePosition).length();

			if (distance < minDistance)
			{
				bestMatch = i;
				minDistance = distance;
				continue;
			}
		}

		toReturn.push_back(bestMatch);
	}
	else if (m_type == EProbeStructureType::Tetrahedral)
	{
	}
	else if (m_type == EProbeStructureType::WeightedNearestNeighbour)
	{
	}

	return toReturn;
}


G3D::Vector3 findNode000(const G3D::Vector3& wsPos, G3D::Vector3& firstProbePosition, float step)
{
    G3D::Vector3 toReturn;

    toReturn.x = floor((wsPos.x - firstProbePosition.x) / step) * step + firstProbePosition.x;
    toReturn.y = floor((wsPos.y - firstProbePosition.y) / step) * step + firstProbePosition.y;
    toReturn.z = floor((wsPos.z - firstProbePosition.z) / step) * step + firstProbePosition.z;

    return toReturn;
}

int findProbeIndex(G3D::Vector3& probePosition, G3D::Vector3& startingPosition, G3D::Vector3& dimensions, float step)
{
    G3D::Vector3 fprobeSpaceCoords = (probePosition - startingPosition) / step;

    G3D::Vector3 probeSpaceCoords = G3D::Vector3(fprobeSpaceCoords);
    int index = probeSpaceCoords.x * dimensions[1] * dimensions[2] +
        probeSpaceCoords.y * dimensions[2] +
        probeSpaceCoords.z;

    return index;
}

ProbeInterpolationRecord ProbeStructure::getInterpolationProbeIndicesAndWeights(const G3D::Vector3& position)
{
	// TODO FIX THIS
	int step = 1;
	ProbeInterpolationRecord record;
	if (m_type == EProbeStructureType::Trilinear)
	{
		//G3D::Array<G3D::Vector3> coords = getInterpolatingProbesCoords(position, step);

		//// Find how close the sample is to probe 000
		//G3D::Vector3 ratios((position[0] - coords[0][0]) / step,
		//					(position[1] - coords[0][1]) / step,
		//					(position[2] - coords[0][2]) / step);

		//std::vector<float> interpolationWeights;// = getInterpolationWeights(ratios);

		//for (float f : interpolationWeights)
		//{
		//	record.weights.push_back(f);
		//}
		// TODO: fix the fact that this doesnt consider the borders.
        G3D::Vector3 firstProbePosition = G3D::Vector3(m_firstProbePosition);
        G3D::Vector3 probe000Pos = findNode000(position, firstProbePosition, step);
        G3D::Vector3 dimensions = G3D::Vector3(m_dimensions[0], m_dimensions[1], m_dimensions[2]);
        int index = findProbeIndex(probe000Pos, firstProbePosition, dimensions, step);
        G3D::Vector3 offsets = G3D::Vector3(m_dimensions[2] * m_dimensions[1], m_dimensions[2], 1);

        G3D::Vector3 w000 = G3D::Vector3(1.0f, 1.0f, 1.0f) - (position - probe000Pos) / step;

		record.weights.push_back(w000.x * w000.y * w000.z);
		record.weights.push_back(w000.x * w000.y * (1.0f - w000.z));
		record.weights.push_back(w000.x * (1.0f - w000.y) * w000.z);
		record.weights.push_back(w000.x *		(1.0f - w000.y) * (1.0f - w000.z));
		record.weights.push_back((1.0f - w000.x) * w000.y * w000.z);
		record.weights.push_back((1.0f - w000.x) * w000.y * (1.0f - w000.z));
		record.weights.push_back((1.0f - w000.x) * (1.0f - w000.y) * w000.z);
		record.weights.push_back((1.0f - w000.x)* (1.0f - w000.y) * (1.0f - w000.z));

        record.probeIndices.push_back(index);
        record.probeIndices.push_back(record.probeIndices[0] + offsets[2]); 							 // p001
        record.probeIndices.push_back(record.probeIndices[0] + offsets[1]); 							 // p010
        record.probeIndices.push_back(record.probeIndices[0] + offsets[1] + offsets[2]);                 // p011
        record.probeIndices.push_back(record.probeIndices[0] + offsets[0]);                              // p100
        record.probeIndices.push_back(record.probeIndices[0] + offsets[0] + offsets[2]);                 // p101
        record.probeIndices.push_back(record.probeIndices[0] + offsets[0] + offsets[1]);                 // p110
        record.probeIndices.push_back(record.probeIndices[0] + offsets[0] + offsets[1] + offsets[2]);    // p111

		int numProbes = m_dimensions[0] * m_dimensions[1] * m_dimensions[2];
		bool shouldNormalize = false;
		for (int i = 0; i < 8; ++i)
		{
			if ((record.probeIndices[i] >= numProbes) || (record.probeIndices[i] < 0))
			{
                record.weights[i] = 0;
				shouldNormalize = true;
			}
		}

		if (shouldNormalize)
		{
			float weightSum = 0;
			for (int i = 0; i < 8; ++i)
			{
				weightSum += record.weights[i];
			}
			for (int i = 0; i < 8; ++i)
			{
                record.weights[i] /= weightSum;
			}
		}
	}
	else if (m_type == EProbeStructureType::Closest)
	{
		
		record.probeIndices = getInterpolatingProbeIndices(position);

		record.weights.push_back(1.0f);
	}
	else if (m_type == EProbeStructureType::Tetrahedral)
	{
		G3D::Array<int> probeIndices;
		G3D::Array<float> weights;
		(App::instance)->tetrahedralInterpolation(position, &probeIndices, &weights);

		record.probeIndices = probeIndices;
		record.weights = weights;
	}

	else if (m_type == EProbeStructureType::WeightedNearestNeighbour)
	{
		G3D::Array<int> probeIndices;
		G3D::Array<float> weights;
		int p = 2;
		float SumOfInverseDistances = 0;
		for (int i = 0 ; i < probeList.size(); ++i)
		{
			Probe* probe = probeList[i];
			G3D::Vector3& probePosition = probe->getPosition();
			float distanceToProbe = (probePosition - position).squaredMagnitude();
			float sumDenominator = 1.f / pow( distanceToProbe, p / 2.f);
			SumOfInverseDistances += sumDenominator;

			probeIndices.append(i);
			float weight = sumDenominator;
			weights.append(weight);
		}

		//float debugsum = 0;
		for (int i = 0; i < weights.size(); ++i)
		{
			weights[i] *= 1.f / SumOfInverseDistances;
			//debugsum += weights[i];
			//debugPrintf("weights[%d]: %f\n", i, weights[i]);
		
		}
		//debugPrintf("sum: %f\n", debugsum);

		record.probeIndices = probeIndices;
		record.weights = weights;

	}

	return record;
}

// Lifted from Stack overflow http://stackoverflow.com/questions/5207550/in-c-is-there-a-way-to-go-to-a-specific-line-in-a-text-file
std::fstream& GotoLine(std::fstream& file, unsigned int num) {
	file.seekg(std::ios::beg);
	for (unsigned int i = 0; i < num - 1; ++i) {
		file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	return file;
}

std::fstream ProbeStructure::probeListFileHandle(bool reading)
{
    int operation = reading ? std::fstream::in : std::fstream::out;
    std::fstream file((probeStructurePath + "/probeList.txt").c_str(), operation);

    if (!file || !file.is_open())
    {
        if (reading)
        {
            debugPrintf("Failed to open probeList file for reading.");

        }
        else
        {
            debugPrintf("Failed to open probeList file for writing.");
        }
    }

    return file;
}

std::fstream ProbeStructure::infoFileHandle(bool reading)
{
	int operation = reading ? std::fstream::in : std::fstream::out;
	std::fstream file((probeStructurePath + "/info.txt").c_str(), operation);

	if (!file || !file.is_open())
	{
		if (reading)
		{
			debugPrintf("Failed to open info file for reading.");

		}
		else
		{
			debugPrintf("Failed to open info file for writing.");
		}
	}

	return file;
}

std::fstream ProbeStructure::probeCoefficientsFileHandle(int i, bool reading)
{
    int operation = reading ? std::fstream::in : std::fstream::out;
    String id = String(std::to_string(i).c_str());
    std::fstream file((probeStructurePath + "/Coefficients/Coefficients_" + id + ".txt").c_str(), operation);

    if (!file || !file.is_open())
    {
        if (reading)
        {
            debugPrintf("Failed to open coeff file for reading.");

        }
        else
        {
            debugPrintf("Failed to open coeff file for writing.");
        }
    }

    return file;
}

std::fstream ProbeStructure::probeCoefficientsGradientsFileHandle(int i, bool reading)
{
    int operation = reading ? std::fstream::in : std::fstream::out;
    String id = String(std::to_string(i).c_str());
    std::fstream file((probeStructurePath + "/Coefficients/CoefficientsGradients_" + id + ".txt").c_str(), operation);

    if (!file || !file.is_open())
    {
        if (reading)
        {
            debugPrintf("Failed to open coeff gradients file for reading.");

        }
        else
        {
            debugPrintf("Failed to open coeff gradients file for writing.");
        }
    }

    return file;
}


void ProbeStructure::updateProbes(bool updateAll)
{
    std::fstream probeListFile = probeListFileHandle(false);

	for (int i = 0; i < probeList.size(); ++i)
	{
		Probe* p = probeList[i];

		// todo: ugh :(
		G3D::Vector3 probePosition;
		if (updateAll)
		{
			probePosition = p->frame.translation;
		}
		else
		{
			probePosition = p->manipulator->frame().translation;
		}

		std::string toWrite = std::to_string(probePosition[0]) + " " + std::to_string(probePosition[1]) + " " + std::to_string(probePosition[2]) + "\n";
		probeListFile << toWrite;
	}

	probeListFile.close();

	bool result = false;

	if (updateAll)
	{
		bool useNew = true;

		// system call to mitsuba
		std::stringstream args;
		args << m_sceneName.c_str() << " " << m_name.c_str() << " Probes";

		runPythonScriptFromDataFiles("onecamera.py", args.str(), true, true);
	}

	for (int i = 0; i < probeList.size(); ++i)
	{
		Probe* p = probeList[i];

		if (updateAll || p->bNeedsUpdate)
		{
			//if (!updateAll)
			//{
			//	// system call to mitsuba
			//	std::stringstream command;
			//	command << "cmd /c \"cd C:\\libraries\\g3d\\samples\\aTest\\data-files\\scripts";
			//	command << " && renderoneprobe.bat " << m_sceneName.c_str() << " " << m_name.c_str() << " " << i << "\"";

			//	result = runCommand(command.str());
			//}

			// reload new probes textures and coefficients
			//if (result == true)
			{
				p->computeCoefficientsFromTexture(false);
				p->bNeedsUpdate = false;
			}

		}
	}
}

bool ProbeStructure::isOutsideSceneBounds(G3D::Vector3 pos, float tolerance)
{
	return App::instance->m_scene->isOOB(pos, tolerance);
}


void ProbeStructure::applyOffsetToProbes(std::vector<float>& displacements)
{
	bool singleAxisMovement = false;
	for (int counter = 0; counter < displacements.size(); counter += 3)
	{
		G3D::Vector3& displacement = G3D::Vector3(displacements[counter * 3 + 0],
												  displacements[counter * 3 + 1],
												  displacements[counter * 3 + 2]);

		if (displacement == G3D::Vector3::zero())
		{
			counter++;
			continue;
		}

		if (singleAxisMovement)
		{
			if ((abs(displacement)[0] > abs(displacement)[1]) && (abs(displacement)[0] > abs(displacement)[2]))
			{
				displacement[1] = 0;
				displacement[2] = 0;
			}
			else if  (abs(displacement)[1] > abs(displacement)[2])
			{
				displacement[0] = 0;
				displacement[2] = 0;
			}
			else
			{
				displacement[0] = 0;
				displacement[1] = 0;
			}
		}

		//float maxDisplacement = 0.05;

		//if (abs(displacement.magnitude()) > maxDisplacement)
		//{
		//	float factor = displacement.magnitude() / maxDisplacement;
		//	displacement /= factor;
		//}

		Probe* probe = probeList[counter];
		G3D::Vector3 newPosition = probe->frame.translation - displacement;
		
		//if (!isOutsideSceneBounds(newPosition, 0.25))
		{
			probe->frame.translation = newPosition;
			probe->bNeedsUpdate = true;
		}


		counter++;
	}
	//offsetFile.close();
}

void ProbeStructure::displaceProbesWithGradient(std::vector<float>& displacements, float maxStepLength)
{
	for (int counter = 0; counter < probeList.size(); counter++)
	{
		G3D::Vector3& displacement = G3D::Vector3(displacements[counter * 3 + 0],
												  displacements[counter * 3 + 1],
												  displacements[counter * 3 + 2]);
        float maxLength = maxStepLength;
        float length = displacement.length();
        if (length > maxLength)
        {
            float factor = sqrtf(length / maxLength);
            displacement /= factor;
        }
		Probe* probe = probeList[counter];

		// Start by sending the probe to its new location
		G3D::Vector3 newPosition = probe->frame.translation - displacement;
		probe->frame.translation = newPosition;
		probe->manipulator->frame().translation = probe->frame.translation;
		probe->bNeedsUpdate = true;

		// Now use the displacement to compute the new coefficients
		for (int k = 0; k < probe->coeffs.size(); ++k)
		{
			G3D::Vector3* coeffs = &(probe->coeffs[k]);
			
			for (int color = 0; color < 3; ++color)
			{
				G3D::Vector3 grad = probe->coeffGradients[k][color];
				
				(*coeffs)[color] += grad.dot(displacement);
			}
			
		}
		//DEAV3(probe->coeffs);
	}
}


TProbeCoefficients ProbeStructure::interpolatedCoefficients(const G3D::Vector3& position, const G3D::Vector3& normal)
{
	ProbeInterpolationRecord record = getInterpolationProbeIndicesAndWeights(position);

	TProbeCoefficients interpolatedCoeffs;

	for (int p = 0; p < record.probeIndices.size(); ++p)
	{
		int index = record.probeIndices[p];
		float weight = record.weights[p];

		TProbeCoefficients& probeCoeffs = probeList[index]->coeffs;
		for (int c = 0; c < m_NumCoefficients; ++c)
		{
			Vector3& coeffs = probeCoeffs[c];

            interpolatedCoeffs.push_back(Vector3());
			interpolatedCoeffs[c].x += weight * coeffs.x;
			interpolatedCoeffs[c].y += weight * coeffs.y;
			interpolatedCoeffs[c].z += weight * coeffs.z;
		}
	}
	

	return interpolatedCoeffs;
}

G3D::Vector3 ProbeStructure::reconstructSH(const G3D::Vector3& position, const G3D::Vector3& normal)
{
    ProbeInterpolationRecord record = getInterpolationProbeIndicesAndWeights(position);
    Vector3 rgb = Vector3(0, 0, 0);

    for (int i = 0; i < record.weights.size(); ++i)
    {
        //int probeIndex = getProbeStartIndex(coords[i], firstProbePosition, dimensions, 1);
        int probeIndex = record.probeIndices[i];
        float weight = record.weights[i];
        if (weight == 0)
        {
            continue;
        }

        //if (this->outputToLog)
        //{
        //    logFile << "Probe " << i << std::endl;
        //    //logFile << "position =  " << coords[i].toString().c_str() << std::endl;
        //    logFile << "weight =  " << weight << std::endl;
        //}

        // For all SH bands
        for (int coeff = 0; coeff < m_NumCoefficients; ++coeff)
        {
            std::pair<int, int> lm = SH::kToLM(coeff);

            // We multiply the inter weight by the geometric term and the SH function value in this direction
            float phong = phongCoeffs(lm.first, 1.0f);
            float sh = SH::SHxyz_yup(lm.first, lm.second, normal);
            float factors = weight *
                //NDotOmegaCoeff(lm(0)) *
                phong *
                sh;
            Vector3& probeCoeffs = getProbe(probeIndex)->coeffs[coeff];
            Vector3 accumulation = factors * probeCoeffs;
            rgb += accumulation;

            //if (this->outputToLog)
            //{
            //    logFile << "SH band: ( " << lm.first << ", " << lm.second << " )" << std::endl;
            //    logFile << "phong = " << phong << std::endl;
            //    logFile << "shFunctionEvaluation = " << sh << std::endl;
            //    logFile << "weight * phong * sh = " << factors << std::endl;
            //    logFile << "probeCoeffs (R,G,B) = " << probeCoeffs.toString().c_str() << std::endl;
            //    logFile << "weight * phong * sh * probeCoeffs = " << accumulation.toString().c_str() << std::endl;
            //    logFile << "accumulatedRGB(" << coeff << ") = " << rgb.toString().c_str() << std::endl << std::endl;
            //}

        }
    }
    rgb /= 3.141592654f;

    return rgb;
}

void ProbeStructure::generateProbes(std::string type)
{
    // system call to mitsuba
    std::stringstream args;
    args << m_sceneName.c_str() << " " << m_name.c_str() << " " << type;

    runPythonScriptFromDataFiles("onecamera.py", args.str(), true, true);

    for (Probe* p : probeList)
    {
        p->texture.reset();
        p->texture = NULL;
    }

}

void ProbeStructure::extractSHCoeffs()
{
    for (int i = 0; i < probeList.size(); ++i)
    {
        Probe* p = probeList[i];

        p->computeCoefficientsFromTexture(true);
		p->saveCoefficients();

		p->bNeedsUpdate = false;
    }

    uploadToGPU();
}

void ProbeStructure::savePositions(bool useManipulator)
{
    std::fstream probeListFile = probeListFileHandle(false);

    for (int i = 0; i < probeList.size(); ++i)
    {
        Probe* p = probeList[i];
        G3D::Vector3 probePosition;
        // todo: ugh :(
        if (useManipulator)
        {
            probePosition = p->manipulator->frame().translation;
            p->setPosition(probePosition);
        }
        else
        {
            probePosition = p->getPosition();
        }

        std::string toWrite = std::to_string(probePosition[0]) + " " + std::to_string(probePosition[1]) + " " + std::to_string(probePosition[2]) + "\n";
        probeListFile << toWrite;

        // check if coeff file exists
        std::fstream coeffFile = probeCoefficientsFileHandle(i, true);
        if (!(coeffFile.is_open()))
        {
            std::fstream newCoeffFile = probeCoefficientsFileHandle(i, false);
            dumpZerosToFile(newCoeffFile, m_NumCoefficients * m_NumColors);
            newCoeffFile.close();

            std::fstream newCoeffGradientFile = probeCoefficientsGradientsFileHandle(i, false);
            dumpZerosToFile(newCoeffGradientFile, m_NumCoefficients * m_NumColors * 3);
            newCoeffGradientFile.close();
        }
        else
        {
            coeffFile.close();
        }

    }

    probeListFile.close();
}

void ProbeStructure::uploadToGPU()
{
    SProbeStructure probeList;

    for (int i = 0; i < probeCount(); ++i)
    {
        SProbe gpuProbe;
        Probe* cpuProbe = getProbe(i);
        for (int j = 0; j < 9; j += 1)
        {
            gpuProbe.coefficients[(3 * j)] = cpuProbe->coeffs[j].x;
            gpuProbe.coefficients[(3 * j) + 1] = cpuProbe->coeffs[j].y;
            gpuProbe.coefficients[(3 * j) + 2] = cpuProbe->coeffs[j].z;

            gpuProbe.gradients[(9 * j) + 0] = cpuProbe->coeffGradients[j][0].x;
            gpuProbe.gradients[(9 * j) + 1] = cpuProbe->coeffGradients[j][0].y;
            gpuProbe.gradients[(9 * j) + 2] = cpuProbe->coeffGradients[j][0].z;

            gpuProbe.gradients[(9 * j) + 3] = cpuProbe->coeffGradients[j][1].x;
            gpuProbe.gradients[(9 * j) + 4] = cpuProbe->coeffGradients[j][1].y;
            gpuProbe.gradients[(9 * j) + 5] = cpuProbe->coeffGradients[j][1].z;

            gpuProbe.gradients[(9 * j) + 6] = cpuProbe->coeffGradients[j][2].x;
            gpuProbe.gradients[(9 * j) + 7] = cpuProbe->coeffGradients[j][2].y;
            gpuProbe.gradients[(9 * j) + 8] = cpuProbe->coeffGradients[j][2].z;

        }
        gpuProbe.position[0] = cpuProbe->position[0];
        gpuProbe.position[1] = cpuProbe->position[1];
        gpuProbe.position[2] = cpuProbe->position[2];
        gpuProbe.position[3] = 0.0f;
        probeList.probes[i] = gpuProbe;
    }

    for (int i = 0; i < m_dimensions.size(); ++i)
    {
        probeList.dimensions[i] = m_dimensions[i];
    }
    probeList.dimensions[3] = 0;;

    probeList.firstProbePosition[0] = m_firstProbePosition[0];
    probeList.firstProbePosition[1] = m_firstProbePosition[1];
    probeList.firstProbePosition[2] = m_firstProbePosition[2];
    probeList.firstProbePosition[3] = 0.0f;

    probeList.step = m_step;

    // Generate the SSBO holding probe information
    GLuint mySSBO = 0;
    glGenBuffers(1, &mySSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, mySSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(probeList), &probeList, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, mySSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

void ProbeStructure::setIntegrator(String& integrator)
{
	this->m_integrator = integrator;

}

void ProbeStructure::setGamma(float gamma)
{
	this->m_gamma = gamma;

}

void ProbeStructure::setType(String& type)
{
	this->m_type = (EProbeStructureType)(int)(typeMap.findIndex(type));
}

void ProbeStructure::saveInfoFile()
{
	std::fstream infoFile = infoFileHandle(false);

	infoFile << "type"			<< " " << std::string(typeMap[m_type].c_str())	<< std::endl;
	infoFile << "dimensions"	<< " " << probeCount()							<< std::endl;
	infoFile << "integrator"	<< " " << std::string(m_integrator.c_str())		<< std::endl;
	infoFile << "gamma"			<< " " << m_gamma								<< std::endl;
	infoFile << "sampleCount"	<< " " << m_NumSamples							<< std::endl;
	infoFile << "width"			<< " " << m_width								<< std::endl;
	infoFile << "height"		<< " " << m_height								<< std::endl;
}

void ProbeStructure::addProbe(G3D::Vector3& position)
{
	Probe* probe = new Probe(probeCount(), probeStructurePath);

	probe->position = position;
	probe->bNeedsUpdate = true;
	probe->m_sphere = Sphere(position, 0.1f);

	probe->frame = CFrame::fromXYZYPRDegrees(position.x, position.y, position.z, 0, 0, 0);

	probe->manipulator = ProbeManipulator::create();
	probe->manipulator->setFrame(probe->frame);
	probe->manipulator->setControlFrame(probe->frame);
	probe->manipulator->p = probe;

	probeList.push_back(probe);
}


bool ProbeStructure::hasProbes()
{
    return probeCount() > 0;
}

void ProbeStructure::saveCoefficients()
{
	for (int i = 0; i < probeList.size(); ++i)
	{
		probeList[i]->saveCoefficients();
	}
}