#include "SceneSample.h"
#include <fstream>

SceneSample::SceneSample(){
}

SceneSample::SceneSample(Vector3 _position){
	this->position = _position;

}

SceneSample::SceneSample(Vector3 _position, Vector3 _normal){
	this->position = _position;
	this->normal = _normal;

}

SceneSample::SceneSample(Vector3 _position, int _selectedModel, Array<int> _probeIndices, Array<float> _weights,Color4 _referenceColor){
	this->position = _position;
	this->selectedModel = _selectedModel;
	this->probeIndices = _probeIndices;
	this->probeWeights = _weights;
	this->referenceColor = _referenceColor;
}

SceneSample::SceneSample(Vector3 _position, int _selectedModel, Color4 _referenceColor){
	this->position = _position;
	this->selectedModel = _selectedModel;

	this->referenceColor = _referenceColor;
}

SceneSample::SceneSample(Vector3 _position, Vector3 _normal, Vector2 _UV, 
						 int _selectedModel, Array<int> _probeIndices, Array<float> _weights, Color4 _referenceColor){
	this->position = _position;
	this->normal = _normal;
	this->UV = _UV;
	this->selectedModel = _selectedModel;
	this->probeIndices = _probeIndices;
	this->probeWeights = _weights;
	this-> referenceColor = _referenceColor;
}

SceneSample::SceneSample(Vector3 _position, Vector3 _normal, Vector2 _IJ, Vector2 _UV, 
						 int _selectedModel, Array<int> _probeIndices, Array<float> _weights, Color4 _referenceColor){
	this->position = _position;
	this->normal = _normal;
	this->UV = _UV;
	this->selectedModel = _selectedModel;
	this->probeIndices = _probeIndices;
	this->probeWeights = _weights;
	this-> referenceColor = _referenceColor;
	this->IJ = _IJ;
}


void SceneSample::writeToFile(String filePath)
{
    std::fstream posFile;
	posFile.open(filePath.c_str(), std::fstream::out | std::fstream::app);
	
	posFile << selectedModel << "\n";

	posFile << this->position.x << " " <<
			   this->position.y << " " <<
			   this->position.z << " " << "\n";
	
	posFile << this->IJ.x << " " <<
			   this->IJ.y << "\n";

	posFile << this->UV.x << " " <<
			   this->UV.y << "\n";

	posFile << this->normal.x << " " <<
			   this->normal.y << " " <<
			   this->normal.z << " " << "\n";

	posFile << this->probeIndices[0] << " " << 
			   this->probeIndices[1] << " " << 
			   this->probeIndices[2] << " " <<  
			   this->probeIndices[3] << "\n";

	posFile << this->probeWeights[0] <<  " " << 
			   this->probeWeights[1] <<  " " << 
			   this->probeWeights[2] <<  " " << 
			   this->probeWeights[3] << "\n";

	posFile << this->referenceColor.r << " " <<
			   this->referenceColor.g << " " <<
			   this->referenceColor.b << "\n";

	posFile << "\n";
	posFile.close();
}

void SceneSample::writeToFile2(int i, int done){
    std::fstream posFile;
	String filename = "../data-files/temp/surfaceSamples/128x128/surfaceSamples_";
	char str[16];
    sprintf(str, "%d", i);
    // sprintf(str2, "%d_%d_%d", i,j,k);
    filename.append(str);
    filename.append(".txt");      
	posFile.open(filename.c_str(), std::fstream::out | std::fstream::app);
	
	posFile << selectedModel << "\n";

	posFile << this->position.x << " " <<
			   this->position.y << " " <<
			   this->position.z << " " << "\n";

	//posFile << this->probeIndices[0] << " " << 
	//		   this->probeIndices[1] << " " << 
	//		   this->probeIndices[2] << " " <<  
	//		   this->probeIndices[3] << "\n";

	//posFile << this->probeWeights[0] <<  " " << 
	//		   this->probeWeights[1] <<  " " << 
	//		   this->probeWeights[2] <<  " " << 
	//		   this->probeWeights[3] << "\n";
	
	posFile << this->referenceColor.r << " " <<
			   this->referenceColor.g << " " <<
			   this->referenceColor.b << "\n";

	if (done == 0) posFile << "1 0 0" << "\n";
	else if (done == 1) posFile << "0 1 0" << "\n";
	else if (done == 2) posFile << "0 0 1" << "\n";
	posFile << "\n";
	posFile.close();
}

std::string SceneSample::toString() const
{
	std::stringstream ss;
	ss << "# Scene sample" << "\n";
	ss << "# Position: \n" << this->position.x << " " << this->position.y << " " << this->position.z  << "\n";
	ss << "# Normal: \n" << this->normal.x << " " << this->normal.y << " " << this->normal.z  << "\n";

	return ss.str();
}