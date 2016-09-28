#include <fstream>
#include "JScene.h"
//#include "G3D/Array.h"
//#include "G3D/stringutils.h"
JScene::JScene()
{
}
JScene::JScene(G3D::String sceneName)
{
	this->m_name = sceneName;

	load();

	if ((m_minBound == G3D::Vector3::zero()) || (m_maxBound == G3D::Vector3::zero()))
	{
		setBounds();
	}
}

void JScene::setBounds()
{
	G3D::Vector3 min = G3D::Vector3(9999.f, 9999.f, 9999.f);
	G3D::Vector3 max = G3D::Vector3(-9999.f, -9999.f, -9999.f);
	for (std::shared_ptr<ArticulatedModel> model : m_models)
	{
		for (ArticulatedModel::Mesh* mesh : model->meshArray())
		{
			for (CPUVertexArray::Vertex& v : mesh->geometry->cpuVertexArray.vertex)
			{
				G3D::Vector3& vertexPosition = v.position;

				for (int axis = 0; axis < 3; ++axis)
				{
					if (vertexPosition[axis] < min[axis])
					{
						min[axis] = vertexPosition[axis];
						continue;
					}
					else if (vertexPosition[axis] > max[axis])
					{
						max[axis] = vertexPosition[axis];
						continue;
					}
				}
			}
		}
	}
	if (m_minBound == G3D::Vector3::zero())
	{
		m_minBound = min;
	}
	if (m_maxBound == G3D::Vector3::zero())
	{
		m_maxBound = max;
	}
	
}

void JScene::load()
{
	m_sceneLights.clear();
    m_meshShapes.clear();
	m_models.clear();
	m_colors.clear();

	std::fstream sceneInfoFile;
	sceneInfoFile.open(("../data-files/Scenes/" + m_name + "/SceneInfo.txt").c_str(), std::ios::in);
	std::string line;
	G3D::Array<G3D::String> tokens;
	while(std::getline(sceneInfoFile, line))
	{
		tokens.clear();
		debugPrintf("read line \" %s \"\n", line.c_str());
		tokens = stringSplit(G3D::String(line.c_str()), ' ');
		if (tokens[0] == "#")
		{
			if (tokens.size() == 1)
			{
				continue;
			}
			if (tokens[1] == "Lights")
			{
				debugPrintf("Setting lights for scene... \n");

				Vector3 position;
				Color3 color;
				while(std::getline(sceneInfoFile, line))
				{
					debugPrintf("read line \" %s \"\n", line.c_str());
					tokens = stringSplit(String(line.c_str()), ' ');
					if ((position != Vector3::zero()) && (tokens[0] == ""))
					{
						addLight(position, color);
						debugPrintf("Added light, position: %s, color: %s \n", position.toString().c_str(), color.toString().c_str());
					}
					else if (tokens[0] == "position")
					{
						position = Point3( std::stof(tokens[1].c_str()), std::stof(tokens[2].c_str()), std::stof(tokens[3].c_str()));
					}
					else if (tokens[0] == "color")
					{
						color = Color3(std::stof(tokens[1].c_str()), std::stof(tokens[2].c_str()), std::stof(tokens[3].c_str()));
					}
					else
					{
						break;
					}
				}

			}
			else if (tokens[1] == "Models")
			{
				debugPrintf("Setting models for scene... \n");
				String filename;
				Color3 color;
				while(std::getline(sceneInfoFile, line))
				{
					debugPrintf("read line \" %s \"\n", line.c_str());
					tokens = stringSplit(String(line.c_str()), ' ');
					if (tokens[0] == "")
					{
						addModel("../data-files/Scenes/" + m_name + "/objs/" + filename, color);
						debugPrintf("Added model, filename: %s, color: %s \n", filename, color.toString().c_str());
					}
					else if (tokens[0] == "filename")
					{
						filename = tokens[1].c_str();
					}
					else if (tokens[0] == "color")
					{
						color = Color3(std::stof(tokens[1].c_str()), std::stof(tokens[2].c_str()), std::stof(tokens[3].c_str()));
					}
					else
					{
						break;
					}
				}
			}
			else if (tokens[1] == "Parameters")
			{
				debugPrintf("Setting parameters for scene... \n");
				while(std::getline(sceneInfoFile, line))
				{
					debugPrintf("read line \" %s \"\n", line.c_str());
					tokens = stringSplit(String(line.c_str()), ' ');
					if (tokens[0] == "scale")
					{
						m_scale = std::stof(tokens[1].c_str());
						debugPrintf("set scale at %f\n", m_scale);
					}
					else if (tokens[0] == "minBound")
					{
						m_minBound = G3D::Vector3(std::stof(tokens[1].c_str()), std::stof(tokens[2].c_str()), std::stof(tokens[3].c_str()));
						debugPrintf("set minbound at %s\n", m_minBound.toString().c_str());
					}
					else if (tokens[0] == "maxBound")
					{
						m_maxBound = G3D::Vector3(std::stof(tokens[1].c_str()), std::stof(tokens[2].c_str()), std::stof(tokens[3].c_str()));
						debugPrintf("set m_maxBound at %s\n", m_maxBound.toString().c_str());
					}
					else
					{
						break;
					}
				}
			}
			else 
			{
			}
		}
	}
	debugPrintf("Finished loading scene! \n");
}

void JScene::addModel(String filename, Color3 color)
{
	if (!loadOBJ(filename, m_scale))
	{
		debugPrintf("Couldn't load file %s", filename);
		return;
	}
	m_colors.push_back(Vector3(color));
}

bool JScene::loadOBJ(String filename)
{
	return loadOBJ(filename, 1.0f);
}

// Add an OBJ file to the sceneModels array and also create its meshshape for sampling
bool JScene::loadOBJ(String filename, float scale)
{
    ArticulatedModel::Specification spec = ArticulatedModel::Specification();
    ParseOBJ parseData;

    spec.filename = System::findDataFile(filename);

	if (spec.filename == "")
	{
		return false;
	}

	spec.scale = scale;
	shared_ptr<ArticulatedModel> model = ArticulatedModel::create(spec);
    m_models.append(model);

	BinaryInput bi(spec.filename, G3D_LITTLE_ENDIAN);
    parseData.parse(bi, spec.objOptions);
    m_meshShapes.append(new MeshShape(parseData));

	return true;
}
void JScene::addLight(Point3 pos, Color3 color)
{
	SceneLight light;
	light.position = pos;
	light.color = color;
	m_sceneLights.append(light);
}

void JScene::addLight(SceneLight light)
{
	m_sceneLights.append(light);
}

int JScene::numModels()
{
	return m_models.size();
}

MeshShape* JScene::getMeshShape(int i)
{
	return m_meshShapes[i];
}

float JScene::scale()
{
	return m_scale;
}

bool JScene::isOOB(G3D::Vector3 pos, float tolerance)
{
	for (int axis = 0; axis < 3; ++axis)
	{
		if ((pos[axis] > (m_maxBound[axis] - tolerance)) || (pos[axis] < (m_minBound[axis] + tolerance)))
		{
			return true;
		}
	}
	return false;
}