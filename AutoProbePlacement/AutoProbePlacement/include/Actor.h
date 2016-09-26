#ifndef ACTOR_H
#define ACTOR_H

#include <G3D/G3DAll.h>

class Actor {
protected:
	ArticulatedModel::Specification spec;
	shared_ptr<ArticulatedModel> model;
	CFrame frame;

	float scale;


public:
  	Array<Vector3> coefficients;
  	float phongExponent;
  	float multiplier;
  	int drawBand;
  	int tetrahedronIndex;
  	bool showOnlyOneBand;
  	bool isVisible;
  	bool updateVisibility;
	shared_ptr<ThirdPersonManipulator> manipulator;
  	String name;
	Actor();
	Actor(String _name, shared_ptr<ArticulatedModel> _model, Point3 translation, shared_ptr<Texture> _texture, bool useManipulator);
	Actor(String name, ArticulatedModel::Specification, Point3, shared_ptr<Texture>);
	shared_ptr<Texture> texture;

	shared_ptr<ArticulatedModel> getModel();
	CFrame getFrame();
	shared_ptr<ThirdPersonManipulator> getManipulator();
	Point3 getPosition();
	void updateCoefficients(Array<Vector3>);
	bool getVisibility();
	void setVisibility(bool b);
	void setPosition(Point3 newPosition);
	Vector3 velocity;
	Vector3 albedo;
};
#endif