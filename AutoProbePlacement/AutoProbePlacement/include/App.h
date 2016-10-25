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
#include "tetgen.h"
#include "ProbeStructure.h"
#include "SceneSample.h"
#include "JScene.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <fstream>

class SceneSampleSet;


struct SProbeOptimization
{
    int id;
    std::vector<float> errors;
    std::vector<float> dp;

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
	G3D::String gamma;
	int filmTypeIndex;
	int integratorIndex;
	bool showWindow;
	bool requireToCloseWindow;
};

struct SProbeStructureCreationOptions
{
	String name;
	String type;
	float step;
	Vector3 firstProbePosition;
	Vector3 dimensions;
};

struct SProbeStructurePanel
{
	int typeIndex;
	int integratorIndex;
    String gamma;
    String width;
    String height;
    String numSamples;
};


enum EFilmType
{
	HDR = 0,
	LDR,
	NUM_FILM_TYPES = LDR + 1
};

class App : public GApp {
private:
	/**
	* METHODS
	*
	* Main app
	*/
  
	void drawModel(RenderDevice* rd, String shaderName, shared_ptr<ArticulatedModel> model, CFrame frame,Args args);

	virtual void makeGui();
    void createRenameOptimizationWindow();
    void createNewProbeWindow();
	void createNewSampleSetWindow();
	void createNewProbe(String& sceneName, String& probeStructureName, G3D::Vector3& position);
	void createNewSampleSet(String& sceneName, String& sNewSampleSetName);
	void createNewProbeStructureWindow();
    void createNewProbeStructure(String& sceneName, String& newSampleSetName);
    void createCopyProbeStructureWindow();
    void copyProbeStructure(String& sceneName, String& sourceProbeStructureName, String& newProbeStructureName);
    void createNewOptimizationSettings();

	ScenePane scenePane;
	GuiPane* probeStructurePane;
	void updateProbeStructurePane();

	void generateSampleSetList();
	void addSampleSetPane(GuiTabPane* tabPane);
	void updateSelectedScenePane();
	void updateProbeStructure();
	void updateSampleSet();
	void addLight(Point3 pos, Color3 color);
	void addLight(SceneLight light);
	void addModel(String filename, Color3 color);
	
	void addActor(String name, String filename, Point3 position, float scale, shared_ptr<Texture> texture);
	void addActor(String name, shared_ptr<ArticulatedModel> model, Point3 position, float scale, shared_ptr<Texture> texture, bool useManipulator, Vector3 albedo);
	
	void debugFloatArray(Array<float> myArray);
	void debugStringArray(Array<String> myArray);
	void debugIntArray(Array<int> myArray);

	void clearAllActors();
	String printVector3Array(Array<Vector3> arr);
	void initializeProbeStructure(String sceneName, String probeStructurePath);
	void loadPreviousProbeStructure();
	void loadScene(String sceneName);

    G3D::String scenesPath();
    G3D::String selectedSceneName();
    G3D::String optimizationFolderPath();
    G3D::String currentOptimizationFolderPath();
    G3D::String probeStructureFoldersPath();
    G3D::String sampleSetFoldersPath();
    G3D::String loadedProbeStructurePath();

    bool probeStructureLoaded();

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


	Vector3 triangleInterpolateVector3(Vector3 weights, Vector3 v0, Vector3 v1, Vector3 v2);
	/**
	* Logic Stuff
	*/

	Array<Vector3> getRandomPoint(int modelNumber, Vector3* P, Vector3* N, Vector3* barycentricWeights, int* startingIndex);
  
	SceneSample generateSceneSample(int* _selectedModel, Vector3 *_P, Vector2* _IJ, Vector3* _N, int sampleID);

	void computeSamplesRGB();
	void computeSamplesRGBRef();
	void computeTriplets();
	void tryOptimization();
	void startOptimizationPasses();
	float computeError();
	void createTempProbeStructure(G3D::Array<G3D::Vector3>& probePositions);
	G3D::Array<G3D::Vector3> generateRandomPositions(int NumberOfPositions);
	void findBestInitialConditions();
	void displaceProbes();

	void addOneActor();
	void addOneActorSq();

	Vector3 getRandomPointInScene();
	bool isUVInsideTriangle(Tri &t, Vector2 testUV, Vector3& barycentric);
	Vector3 getUVBarycentricCoordinates(Point2 testUV, Vector3 u, Vector3 v);
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


	/**
	* Various
	*/
	// Number of bands we will use for SH reconstruction
	int numBands;
	int totalNumberOfCoeffs;
	int maxDrawBand;
	int drawBand;

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
    

	std::string tetgenString;

	Array<Vector3> randomPoints;
	Vector3* rndPtBaryWeight;
	int* startingIndex;
	int selectedModel;
	String *samplesToSave;
	String* SampleSetOutputName;

	Random m_random;

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
	bool bOptimizeForCoeffs;
	TriTree triTree;
	//Index* trisIndex;

public:
	float shadingMultiplier;
	ProbeStructure *m_probeStructure;

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
    virtual void onUserInput(UserInput* userInput);
	String m_scenePath;
	G3D::String previousProbeStructure;
	G3D::String tbNumPassesLeft;
    G3D::String numOptimizationSamples;

	int numPassesLeft;
	bool logSampleSet;

	Array<Vector3> extrapolationTriangleVertices;
	void tetrahedralInterpolation(G3D::Vector3 testPoint, Array<int> *_probeIndices, Array<float> *_coeffs);

	void switchEditProbeStructure();
	JScene* m_scene;

	String m_sNumICTries;
	String m_sNumICProbes;

	void offlineRender();

	void clearSampleSetValues();
	void clearSampleSetPositions();
	void reloadSampleSet();
	void generateSampleSetPositions();
	void generateSampleSetValues();
    void generateSampleSetValuesFromProbes();

	// There CERTAINLY has to be a better way to handle this stuff than keeping the new sample set window and nmae as globals
	G3D::String sNewSampleSetName;
	G3D::String sNewProbePosition;
	SProbeStructureCreationOptions newProbeStructureOptions;


	shared_ptr<GuiWindow> windowNewProbe;
	shared_ptr<GuiWindow> windowNewSampleSet;
    shared_ptr<GuiWindow> windowNewProbeStructure;
    shared_ptr<GuiWindow> windowCopyProbeStructure;
    shared_ptr<GuiWindow> windowRenameOptimization;
	SOfflineRenderingOptions offlineRenderingOptions;
	SProbeStructurePanel probeStructurePanelOptions;

	G3D::Array<G3D::String> integratorList;
	G3D::Array<G3D::String> filmTypeList;

    int m_currentOptimization;
    String sCurrentOptimizationName;
    SProbeOptimization currentOptimization;

    bool bShouldUpdateProbeStructurePane;
    bool bKeepRefValuesOnNewOptimization; 
    bool bTakeRefScreenshot;
    bool bOptimizeWithMitsubaSamples;
};

#endif