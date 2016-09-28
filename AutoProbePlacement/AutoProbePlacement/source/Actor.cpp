#include "Actor.h"

Actor::Actor(){

}
Actor::Actor(String _name, ArticulatedModel::Specification _spec, Point3 translation, shared_ptr<Texture> _texture){
	//
	this->spec = _spec;
	this->model = ArticulatedModel::create(_spec);

	this->frame = CFrame::fromXYZYPRDegrees(translation.x,translation.y,translation.z,0,0,0);


	this->manipulator = ThirdPersonManipulator::create();
	this->manipulator->setFrame(frame);
    this->manipulator->setControlFrame(frame);

    this->multiplier = 1.0;
    this->phongExponent = 1.0;
    this->drawBand = 8;
    this->showOnlyOneBand = false;
    this->tetrahedronIndex = 5;

    this->isVisible = true;
    this -> updateVisibility = false;

    this->name = _name;

    this->texture = _texture;

    float x,y,z;
    Random::common().sphere(x,y,z);
    this->velocity = Vector3(x,y,z);

	for (int i = 0; i < 9; ++i)
	{
		this->coefficients.append(Vector3(0,0,0));
	}
}

Actor::Actor(String _name, shared_ptr<ArticulatedModel> _model, Point3 translation, shared_ptr<Texture> _texture, bool useManipulator){

	this->model = _model;

	this->frame = CFrame::fromXYZYPRDegrees(translation.x,translation.y,translation.z,0,0,0);


	if (useManipulator)
	{
		this->manipulator = ThirdPersonManipulator::create();
		this->manipulator->setFrame(frame);
		this->manipulator->setControlFrame(frame);
	}
	else
	{
		this->manipulator = NULL;
	}

    this->multiplier = 1.0;
    this->phongExponent = 1.0;
    this->drawBand = 8;
    this->showOnlyOneBand = false;
    this->tetrahedronIndex = 0;

    this->isVisible = true;
    this -> updateVisibility = false;

    this->name = _name;

    this->texture = _texture;

    float x,y,z;
    Random::common().sphere(x,y,z);
    this->velocity = Vector3(x,y,z);

	for (int i = 0; i < 9; ++i)
	{
		this->coefficients.append(Vector3(0,0,0));
	}
}

shared_ptr<ArticulatedModel> Actor::getModel(){
	return model;
}

CFrame Actor::getFrame(){
	return this->manipulator->frame();
}

shared_ptr<ThirdPersonManipulator> Actor::getManipulator(){
	return manipulator;
}

void Actor::updateCoefficients(Array<Vector3> _coefficients){
	this->coefficients = _coefficients;
}

Point3 Actor::getPosition(){
	return this->manipulator->frame().translation;
}

void Actor::setPosition(Point3 newPosition){
	// logPrintf("trying to set actor position to %s. before: %s\n",newPosition.toString().c_str(), this->manipulator->frame().toXYZYPRDegreesString().c_str());
	this->frame = CFrame(this->manipulator->frame().rotation, newPosition);
	this->manipulator->setFrame(frame);
    // this->manipulator->setControlFrame(frame);
	// logPrintf("after: %s\n", this->manipulator->frame().toXYZYPRDegreesString().c_str());
}

bool Actor::getVisibility(){
	return this->isVisible;
}
void Actor::setVisibility(bool b){
	updateVisibility = b;
}

