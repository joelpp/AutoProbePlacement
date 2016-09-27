/**
@file App.h

The G3D 8.0 default starter app is configured for OpenGL 3.0 and relatively recent
GPUs.  To support older GPUs you may want to disable the framebuffer and film
classes and use G3D::Sky to handle the skybox.
*/
#ifndef App_h
#define App_h

#include <G3D/G3DAll.h>
#include <GLG3D/GLG3D.h>
#include "Probe.h"
#include "Actor.h"
#include "Tetrahedron.h"
#include "Transform.h"
#include "spline.h"
#include "tetgen.h"
#include "ProbeStructure.h"
#include "SceneSample.h"
#include "JScene.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <fstream>

class SceneSampleSet;

struct SProbe
{
	float position[4];
	float coefficients[27];
	float gradients[81];
};

struct SProbeStructure
{
	SProbe probes[343];
	int dimensions[4]; // int3 padded into int4
	float firstProbePosition[4];  // vec3 padded into vec4
	float step;
};


struct ScenePane
{
	G3D::GuiDropDownList* selectedSceneList;
	G3D::GuiDropDownList* probeStructureList;
	G3D::GuiDropDownList* sampleSetList;
};

struct SOfflineRenderingOptions
{
	G3D::String numSamples;
	G3D::String width;
	G3D::String height;
	bool showWindow;
	bool requireToCloseWindow;
};



class App : public GApp {
private:
	SProbeStructure probeList;
	/**
	* METHODS
	*
	* Main app
	*/
  
	void drawModel(RenderDevice* rd, String shaderName, shared_ptr<ArticulatedModel> model, CFrame frame,Args args);
	void makeScene();
	virtual void makeGui();
	void addScenePane(GuiTabPane* tabPane);
	void updateSelectedScenePane();
	void updateProbeStructure();
	void updateSampleSet();
	void addLight(Point3 pos, Color3 color);
	void addLight(SceneLight light);
	void addModel(String filename, Color3 color);
	
	void makeZTextures();
	void loadSponza();
	void loadZScene();
	void loadSimpleTriangle();
	void makeBoxTextures();
	void makeSponzaTextures();
	void loadCornell();
	void loadBoxScene();
    void loadCornellColors();

	void loadSHTextures();
	void saveCoeffs(String path);
	void addActor(String name, String filename, Point3 position, float scale, shared_ptr<Texture> texture);
	void addActor(String name, shared_ptr<ArticulatedModel> model, Point3 position, float scale, shared_ptr<Texture> texture, bool useManipulator, Vector3 albedo);
	void addParticle(float x, float y, float z);
	void addParticle();
	void resetParticles();
	void debugFloatArray(Array<float> myArray);
	void debugStringArray(Array<String> myArray);
	void debugIntArray(Array<int> myArray);
	//bool loadOBJ(String filename);
	//bool loadOBJ(String filename, float scale);
	void clearAllActors();
	void pythonTest();
	void updateSelectedSceneTextures();
	String printVector3Array(Array<Vector3> arr);
	void logPrintSurfaceSamples();
	void initializeProbeStructure(String sceneName, String probeStructurePath);
	void loadPreviousProbeStructure();
	void loadScene(String sceneName);

	/**
	* Interpolation
	*/
	Array<Vector3> getInterpolatedCoeffs(ProbeInterpolationRecord iRec, int maxBand);
	Array<Vector3> getTriangleVertices(double T, Tetrahedron* tet);

	Vector3 getBarycentricCoordinates(Point3 testPoint, Tetrahedron* t);
	Vector3 findTriangleBarycentricCoordinates(Vector3 P, Vector3 T0, Vector3 T1, Vector3 T2);

	Tetrahedron* updateNeighbor(Vector3 weights, float d, Actor& actor, Tetrahedron* tet);

	bool tetrahedralInterpolation(Actor& actor, Array<int> *_probeIndices, Array<float> *_coeffs);
	bool badWeights(Vector3 v, float d);

	double findCubicRoot(double a0, double a1, double a2, double a3);
	double findTforCurrentTriangleFace(Point3 testPosition, Tetrahedron* tet);
	double extrapolationMatrixDeterminantFunction(double t, Vector3 V0, Vector3 V1, Vector3 V2, Vector3 N0, Vector3 N1, Vector3 N2);
	double extrapolationMatrixDeterminantFunctionDerived(double t, Vector3 V0, Vector3 V1, Vector3 V2, Vector3 N0, Vector3 N1, Vector3 N2);

	int minCoord(Vector3 weights);

	void setToCorrectTetrahedron(int actorID);

	Array<Vector3> triLinearInterpolation(int actorID, Array<int> *_probeIndices, Array<float> *_coeffs, bool keepProbes);
	Vector3 triangleInterpolateVector3(Vector3 weights, Vector3 v0, Vector3 v1, Vector3 v2);
	/**
	* Logic Stuff
	*/

	Array<Vector3> getRandomPoint(int modelNumber, Vector3* P, Vector3* N, Vector3* barycentricWeights, int* startingIndex);
  
	Vector3 SHDotProductWithCosineLobe(Array<Vector3> coefficients);
	Vector3 shCoeffsToRGB(Vector3 normal, Array<Vector3> coeffs);

	void sampleScenePoint(int* _selectedModel, Vector3 *_P, Vector2* _IJ, Vector3* _N, int sampleID);

	void computeSceneSamples();
	void computeSamplesRGB();
	void computeSamplesRGBRef();
	void computeTriplets();
	void tryOptimization();
	void startOptimizationPasses();
	float computeError(std::string logFilePath);
	void createTempProbeStructure(G3D::Array<G3D::Vector3>& probePositions);
	G3D::Array<G3D::Vector3> generateRandomPositions(int NumberOfPositions);
	void findBestInitialConditions();
	void displaceProbes();

	void reloadProbeStructure();
	void logicFunWithCoeffs();
	void logicFunWithPositions();
	void callLogicFunNTimes();
	void addOneActor();
	void addOneActorSq();
	void extractProbeCoeffs();
	void probePositionTesting();
	void probeStructureTesting();
	void callLogicFunNTimesPos();

	float clampedCosLobeCoefficient(int l);
	float phongCoeffs(int l, float r);

	Array<int> kToLM(int k);
	bool atLeast2DifferentCoordinates(Point3 p0, Point3 p1);
	Array<Vector3> findOppositeCornersOnRectangle(Array<Point3> vertexCoords);
	Array<Vector3> getModelCorners(int selectedModel);
	Array<Vector3> findModelOppositeCorners(int i);
	void discretizeModel(int selectedModel);
	Vector3 getRandomPointInScene();
	void loadSurfaceSamples(String& sceneName, String& sampleSetName);
	void loadSampleCoeffs();
	Point3 findMinimalestCorner(Array<Point3> vertexCoords);
	bool isUVInsideTriangle(Tri &t, Vector2 testUV, Vector3& barycentric);
	Vector3 getUVBarycentricCoordinates(Point2 testUV, Vector3 u, Vector3 v);
	void discretizeScene();
	bool isUVInsideModel(Vector2 UV, Tri &tri, Vector3 &barycentric);
	void testTriangleFunctions();
	bool pointInTriangleBB(Point2 UV, Vector3 u, Vector3 v);

	////////////////////////////////
	///////// RENDERING ////////////
	////////////////////////////////

	void renderActors(RenderDevice* rd);
	void drawProbes(RenderDevice* rd);
	void drawLights(RenderDevice* rd);
	void drawScene(RenderDevice* rd);
	void drawSurfaceSamples(RenderDevice* rd);
	void drawProbeLineSegments(RenderDevice* rd);
	void setSHTexturesUniforms(Args& args);
	void setProbeCoeffUniforms(Args& args, G3D::Array<G3D::Vector3>& coeffs);


	/*
	* GLOBAL VARIABLES
	*
	* 
	* Textures
	*/
	Array<shared_ptr<Texture> >         shTextures;
	Array<shared_ptr<Texture> >         bakedSceneTextures;
	Array<shared_ptr<Texture> >         optimizedSceneTextures;

	Array<Array<shared_ptr<Texture> > > sceneTextures;

	shared_ptr<Texture> billiardTexture;

	/**
	* Models
	*/
	shared_ptr<ArticulatedModel>        dragonModel;
	shared_ptr<ArticulatedModel>        sceneModel;
	shared_ptr<ArticulatedModel>        sphereModel;
	shared_ptr<ArticulatedModel>        squareModel;
	shared_ptr<ArticulatedModel>        bunnyModel;
	shared_ptr<ArticulatedModel>        model;
	shared_ptr<ArticulatedModel>        currentModel;
	shared_ptr<ArticulatedModel> armadilloModel;
	CPUVertexArray cva, cvaUV;
	Array<Tri> triArray, triArrayUV;;
	/**
	* Rendering
	*/
	shared_ptr<ThirdPersonManipulator>  manipulator;
	Array<shared_ptr<Surface> >         m_sceneGeometry;
	LightingEnvironment                 environment;
	AttributeArray                      attArray;
	Array<AttributeArray> shAttributeArray;
	Array<AttributeArray> colorsAtt;
	Array<Vector3> coeffValues;

	ProbeStructure *OGprobeStructure;
	ProbeStructure *anotherProbeStructure;
	Array<Probe*> probesToRender;
	CFrame probeTransform;
	Actor particle;
	Vector3 velocity;
	int timer;
	Array<Actor> actors;
	Array<Actor> particles;
	float lightIntensity;
	float sceneTextureIntensity;

	Array<LineSegment> lineArray;
	Array<Tetrahedron*> tetrahedraArray;
	Array<Array<int> > tetFaces;
	Array<Vector3> tetFaceNormals;
	Array<Vector3> tetFaceNormalPositions;
	double extrapolationT;
	int smallestBaryCoord;

	bool useIrradianceTexture;
	bool bManipulateProbesEnabled;

	ScenePane scenePane;

	/**
	* Various
	*/
	// Number of bands we will use for SH reconstruction
	int numBands;
	int totalNumberOfCoeffs;
	int maxDrawBand;
	int drawBand;

	float shadingMultiplier;
	float sampleMultiplier;
	float phongExponent;
	float yaw;
	float pitch; 
	float roll;
  
	Point3 spherePosition;
	Point3 lightPosition;
	Array<Point3> lightPositions;
	Array<SceneLight> sceneLights;

	bool showInterpolationProbes, showAllProbes;
	bool CPUInterpolation;
	bool highlightProbes;
	bool init;
	bool showOnlyOneBand;
	bool showFaceNormals;
	bool showParticles;
	bool hideCeiling;
	bool hideRoom;
	bool hideActors;
	bool highlightTetrahedron;
	bool saveSample;
	bool useBakedSceneTextures;
	bool interpolateCoefficients;


	Vector3 offset;
	Vector3 step;

	std::string tetgenString;

	Array<Vector3> randomPoints;
	Vector3* rndPtBaryWeight;
	int* startingIndex;
	int selectedModel;
	String *samplesToSave;
	String* SampleSetOutputName;

	Random r;

	float prevCost;

	String actorSpawnX, actorSpawnY, actorSpawnZ;

	String OSXPATH;

	Array<Point3> modelPoints;
	Array<Point3> samplePoints;
	Array<SceneSample> sceneSamples;
	SceneSampleSet* sampleSet;
	int sampleDropDownIndex;
	String maxSamplesPointsToDraw;
	bool showSampleNormals;
	bool showDarkSamples;
	bool showSamples;
	bool useSHGradients;
	bool useMatlabOptimization;
	TriTree triTree;
	//Index* trisIndex;

public:
	ProbeStructure *probeStructure;

	bool bRenderDirect;
	bool bRenderIndirect;
	bool bRenderMultiplyIndirectByBRDF;
	bool bRenderShadowMaps;
	bool bRenderAO;

	static App* instance;
	App(const GApp::Settings& settings = GApp::Settings());
	virtual void onAI() override;
	virtual void onInit();
	virtual void onGraphics3D(RenderDevice* rd, Array< shared_ptr<Surface> >& surface);
	String m_scenePath;
	G3D::String previousProbeStructure;

	int numPassesLeft;
	G3D::String tbNumPassesLeft;
	bool logSampleSet;

	Array<Vector3> extrapolationTriangleVertices;
	void tetrahedralInterpolation(G3D::Vector3 testPoint, Array<int> *_probeIndices, Array<float> *_coeffs);

	void switchEditProbeStructure();
	void saveProbeStructureUpdateAll();
	void saveProbeStructureUpdate();
	JScene m_scene;

	String m_sNumICTries;
	String m_sNumICProbes;

	void offlineRender();

};

#endif