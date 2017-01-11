#ifndef SCENE_H
#define SCENE_H

//#include <memory>
#include <string>
#include <G3D/G3DAll.h>

struct SceneLight
{
	G3D::Point3 position;
	G3D::Color3 color;
};
class JScene {
protected:


public:

	JScene();
	JScene(G3D::String sceneName);
	void load();
	void addModel(String filename, Color3 color);
	bool loadOBJ(String filename);
	bool loadOBJ(String filename, float scale);
	void addLight(Point3 pos, Color3 color);
	void addLight(SceneLight light);
	void setBounds();
	int numModels();
	MeshShape* JScene::getMeshShape(int i);
	float scale();
	bool isOOB(G3D::Vector3 pos, float tolerance);

    G3D::String name() { return m_name; }

	G3D::String m_name;
	float m_gradientDisplacement;
	G3D::Vector3 m_minBound;
	G3D::Vector3 m_maxBound;
	float m_scale;
	G3D::Array<G3D::Vector3> m_colors;
	G3D::Array<std::shared_ptr<G3D::ArticulatedModel> > m_models;
	G3D::Array<G3D::MeshShape*> m_meshShapes;
	G3D::Array<SceneLight> m_sceneLights;

};

#endif // SCENE_H