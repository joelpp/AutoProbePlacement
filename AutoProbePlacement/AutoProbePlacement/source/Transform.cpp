#include "Transform.h"

Vector2 Transform::cartesianToSpherical(Vector3 xyz){
	Vector2 PT = Vector2();

	PT.x = atan2(xyz.z, xyz.x);
	PT.y = G3D::acos(xyz.y);

	return PT;
}

Vector3 Transform::sphericalToCartesian(Vector2 PT, float r = 1.0f)
{
	Vector3 xyz;

	float phi = PT.x;
	float theta = PT.y;

	xyz.x = r * sin(theta) * cos(phi);
	xyz.y = r * cos(theta);
	xyz.z = r * sin(theta) * sin(phi);

	return xyz;
}

Vector2 Transform::cartesianToUV(Vector3 xyz)
{
	return sphericalToUV(cartesianToSpherical(xyz));
}

Vector2 Transform::sphericalToUV(Vector2 PT)
{
	float pi = 3.141592654;
	return Vector2(PT.x / pi, 2.0f*PT.y/pi - 1.0f);
}



Vector2 Transform::uvToIJ(Vector2 UV, int width, int height)
{
    Vector2 toReturn = Vector2();

    toReturn.y = height * (UV.y + 1.0) / 2.0;
    toReturn.x = width * (UV.x + 1.0) / 2.0;

    return toReturn;
}

Vector2 Transform::ijToUV(Vector2 IJ, int width, int height)
{
	return Vector2(2.0 * IJ.x / width - 1,  2.0 * IJ.y / height - 1);
}

Vector2 Transform::UVtoPT(Vector2 UV)
{
	float pi = 3.141592654;
	return Vector2(pi*UV.x, pi*(UV.y + 1) / 2.0);
}

Vector3 Transform::PTtoXYZ(float r, Vector2 PT)
{
	return Vector3(r * sin(PT.y) * cos(PT.x) , r * cos(PT.y) ,  r * sin(PT.y) * sin(PT.x));
}
Vector2 Transform::worldSpaceToUV(Vector3 weights, shared_ptr<ArticulatedModel> model, Array<Vector3> vertices)
{
    
   
    Vector2 UV = Vector2();

    Vector2 v0, v1, v2;

    Array<Vector2> texCoords = Array<Vector2>();

    texCoords.append(Vector2());
    texCoords.append(Vector2());
    texCoords.append(Vector2());
    for (int i = 0; i < vertices.size(); i++){
        for (int j = 0; j < model->geometryArray()[0]->cpuVertexArray.size(); j++)
            if (model->geometryArray()[0]->cpuVertexArray.vertex[j].position == vertices[i])
                texCoords[i] = model->geometryArray()[0]->cpuVertexArray.vertex[j].texCoord0;
    }

    UV = texCoords[0] * weights.x + texCoords[1] * weights.y + texCoords[2] * weights.z;


    return UV;
}

Color3 Transform::color3FromVector3(Vector3 toTransform){ return Color3(toTransform.x, toTransform.y, toTransform.z); }
Vector3 Transform::Vector3FromColor3(Color3 toTransform){ return Vector3(toTransform.r, toTransform.g, toTransform.b); }

Vector2 Transform::V2fromString(String s){
	Vector2 toReturn = Vector2();
	Array<String> stringArray = stringSplit(s, ' ');

	return Vector2(0, 0);
}

Vector3 Transform::Vector3FromString(String toTransform){
	
	Array<String> stringArray = stringSplit(toTransform, ' ');

	return Vector3(atof(stringArray[0].c_str()), atof(stringArray[1].c_str()), atof(stringArray[2].c_str()));
}
Color3 Transform::Color3FromString(String toTransform){
	
	Array<String> stringArray = stringSplit(toTransform, ' ');

	return Color3(atof(stringArray[0].c_str()), atof(stringArray[1].c_str()), atof(stringArray[2].c_str()));
}
