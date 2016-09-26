#include <G3D/G3DAll.h>

class Transform {
protected:

public:
	static Vector2 cartesianToSpherical(Vector3 xyz);
	static Vector3 sphericalToCartesian(Vector2 PT, float r);
	static Vector2 cartesianToUV(Vector3 xyz);
	static Vector3 PTtoXYZ(float r, Vector2 PT);
	static Vector2 uvToIJ(Vector2 UV, int width, int height);
	static Vector2 ijToUV(Vector2 IJ, int width, int height);
	static Vector2 UVtoPT(Vector2 UV);
	static Vector2 sphericalToUV(Vector2 PT);
	static Color3 color3FromVector3(Vector3 toTransform);
	static Vector3 Vector3FromColor3(Color3 toTransform);
	static Vector3 Vector3FromString(String toTransform);
	static Color3 Color3FromString(String toTransform);

	static Vector2 worldSpaceToUV(Vector3 weights, shared_ptr<ArticulatedModel> model, Array<Vector3> vertices);
	static Vector2 V2fromString(String s);

};