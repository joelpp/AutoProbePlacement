#ifndef SCENESAMPLE_H
#define SCENESAMPLE_H

#include <G3D/G3DAll.h>
#include <sstream> 
class SceneSample{
	public:
		Vector3 position;
		Vector3 normal;
		Array<float> values;
		Vector2 UV;
		Color4 referenceColor;
		Vector2 IJ;
		int selectedModel;
		Array<int> probeIndices;
		Array<float> probeWeights;
		
		SceneSample(Vector3 _position, int _selectedModel, Color4 _referenceColor);
		SceneSample();
		SceneSample(Vector3 _position, int _selectedModel, Array<int> _probeIndices, Array<float> _weights, Color4 _referenceColor);
		SceneSample(Vector3 _position, Vector3 _normal, Vector2 _UV, int _selectedModel, Array<int> _probeIndices, Array<float> _weights, Color4 uvColor);
		SceneSample(Vector3 _position, Vector3 _normal, Vector2 _IJ, Vector2 _UV, int _selectedModel, Array<int> _probeIndices, Array<float> _weights, Color4 uvColor);
		SceneSample(Vector3 _position, Vector3 _normal);
		SceneSample(Vector3 _position);
		void writeToFile(String filePath);
		void writeToFile2(int i, int done);
		std::string toString() const;
};


#endif