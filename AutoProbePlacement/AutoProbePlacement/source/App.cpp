/** \file App.cpp */
#include <fstream>

#include "App.h"
#include "SceneSampleSet.h"
#include "ProbeRenderer.h"
#include "Helpers.h"

/*
    Debugging Helpers
 */
#define DE_(x) debugPrintf(#x);debugPrintf("\n\n");
#define DES(x) debugPrintf(#x);debugPrintf(": %s\n",x.c_str());
#define DEB(x) debugPrintf(#x);debugPrintf(": %f\n",x)
#define DEB2(x,y) debugPrintf(#x);debugPrintf(": %f\t",x);debugPrintf(#y);debugPrintf(": %f\n",y)
#define DEB3(x,y,z) debugPrintf(#x);debugPrintf(": %f\t",x);debugPrintf(#y);debugPrintf(": %f\t",y);debugPrintf(#z);debugPrintf(": %f\n",z)
#define DEI(x) debugPrintf(#x);debugPrintf(": %d\n",x)
#define DEI2(x,y) debugPrintf(#x);debugPrintf(": %d\t",x);debugPrintf(#y);debugPrintf(": %d\n",y)
#define DEI3(x,y,z) debugPrintf(#x);debugPrintf(": %d\t",x);debugPrintf(#y);debugPrintf(": %d\t",y);debugPrintf(#z);debugPrintf(": %d\n",z)
#define DEV2(v) debugPrintf(#v);debugPrintf(": (%f,\t%f)\n",v.x,v.y)
#define DEV3(v) debugPrintf(#v);debugPrintf(": (%f,\t%f,\t%f)\n",v.x,v.y,v.z)
#define DEV4(v) debugPrintf(#v);debugPrintf(": (%f,\t%f,\t%f,\t%f)\n",v.x,v.y,v.z,v.w)
#define DEC(v) debugPrintf(#v);debugPrintf(": (%f,\t%f,\t%f)\n",v.r,v.g,v.b)
#define DEAV3(x) debugPrintf(#x); for(int num = 0; num < x.size(); num++){ debugPrintf(", [num]: (%s)\n",x[num].toString().c_str()); debugPrintf("\n");}

G3D_START_AT_MAIN();

// Function Declarations
float legendreP(int l_, int m_, float z);
float SH(int m, int l, Point3 pos);
float S(int m, float x, float y);
float C(int m, float x, float y);

float findDistance(Vector3 v1, Vector3 v2);
Vector3 transform(Vector3 old);



int main(int argc, char** argv) {


    (void)argc; (void)argv;
    GApp::Settings settings;
    
    settings.window.width       = 1440; 
    settings.window.height      = 900;
    settings.window.caption     = "Assisted probe placement";
	settings.dataDir = "C:\\git\\AutoProbePlacement\\data-files";
	settings.renderer.deferredShading = true;
	settings.renderer.factory = &(ProbeRenderer::create);
	//settings.film.preferredColorFormats.insert(0,ImageFormat::RGBA32F());
//	settings.screenshotDirectory = "C:/temp/screenshots";
//#   ifdef G3D_WINDOWS
//        // On Unix operating systems, icompile automatically copies data files.  
//        // On Windows, we just run from the data directory.
//        if (FileSystem::exists("data-files")) {
//            chdir("data-files");
//        } else if (FileSystem::exists("../samples/atest/data-files")) {
//            chdir("../samples/atest/data-files");
//        }
//
//#   endif
    return App(settings).run();
}


App::App(const GApp::Settings& settings) : GApp(settings) {}

void App::onInit() {
	GApp::onInit();

	instance = this;
    renderDevice->setSwapBuffersAutomatically(true);
    
    actors = Array<Actor>();

    //Load sphere model    
    ArticulatedModel::Specification spec;
    spec = ArticulatedModel::Specification();
    spec.scale = 1.0f;
	logPrintf("lol");

	// this should be some kind of dict or sane management system
    spec.filename = System::findDataFile("../data-files/objs/sphere.obj");
    spec.stripMaterials = false;
    spec.scale          = 0.1f;
    sphereModel = ArticulatedModel::create(spec);
	spec.filename = System::findDataFile("../data-files/objs/square.obj");
	squareModel = ArticulatedModel::create(spec);
	spec.filename = System::findDataFile("../data-files/objs/bunny.obj");
	bunnyModel = ArticulatedModel::create(spec);

    spec.scale = 1.0f;
	//spec.filename = System::findDataFile("../data-files/Scenes/crytek-sponza/sponza_v2_scaled.obj");
	//spec.filename = System::findDataFile("../data-files/Scenes/ZScene/objs/largegreenwall.obj");
	//sceneModel = ArticulatedModel::create(spec);


    //Set Critical settings
    // spherePosition = Point3(0,4.5,0);
    lightPosition = Point3(0,4,0);
    // lightPosition = Point3(0,2.5,0);
    showInterpolationProbes = false;
    showAllProbes = false;
    CPUInterpolation = false;
    highlightProbes = false;
    showFaceNormals = false;
    hideCeiling = true;
    hideRoom = false;
    showParticles = false;
	useIrradianceTexture = false;
	highlightTetrahedron = false;
	saveSample = false;
	useBakedSceneTextures = true;
	interpolateCoefficients = false;
	logSampleSet = false;
    lightIntensity = 1;
    sceneTextureIntensity = 1.0;
	maxSamplesPointsToDraw = "0";
	showSampleNormals = false;
	showDarkSamples = true;
	showSamples = true;
	useSHGradients = false;
	bManipulateProbesEnabled = false;
	hideActors = false;
	useMatlabOptimization = true;
	bRenderDirect = true;
	bRenderIndirect = true;
	bRenderMultiplyIndirectByBRDF = false;
	sampleSet = NULL;
    timer = 0;
	m_activeCamera->setPosition(Point3(0,3,3));
	sampleDropDownIndex = 0;
	tetgenString = "../data-files/tetgen/44/";
	previousProbeStructure = "";

	samplesToSave = new String("1");
	SampleSetOutputName = new String("SampleSetOutputName");
    //Decide how many bands we want to use for the interpolation
    numBands = 2;
    totalNumberOfCoeffs = 0;
    maxDrawBand = 8;
	shadingMultiplier = 1.0f;
	sampleMultiplier = 1.0f;
    for (int i = 0; i <= numBands; i++)
	{
		totalNumberOfCoeffs += 2*i + 1; 
	}

    extrapolationT = 0;
    smallestBaryCoord = -1;
	actorSpawnX = String("0");
	actorSpawnY = String("0"); 
	actorSpawnZ = String("0");
	tbNumPassesLeft = String("1");
	numPassesLeft = 0;
    createDeveloperHUD();

    //Add sliders for rendering parameters
	addOneActor();
    makeGui();

	String sceneName;
	String sampleSetName;


	probeStructure = NULL;

    //Load the textures containing the SH coefficient values
    //loadSHTextures();
	GApp::loadScene((String)"C:\\git\\AutoProbePlacement\\AutoProbePlacement\\data-files\\Scenes\\zcbox\\cbox.Scene.Any");

}//end of onInit 

void App::loadScene(String sceneName)
{
	lightPositions.clear();

	if (sceneName == "ZScene")
	{
		loadZScene();
	}
	else if (sceneName == "crytek-sponza")
	{
		loadSponza();
	}
	//else if (sceneName == "zcbox")
	//{
	//	loadCornell();
	//}
	else
	{
		m_scene = JScene(sceneName);
		//loadScene(std::string(sceneName.c_str()));
	}

}

void App::initializeProbeStructure(String sceneName, String probeStructureName)
{
	//Load the probes from the currently active experiment
	String probeStructurePath = "../data-files/Scenes/" + sceneName + "/ProbeStructures/" + probeStructureName;

	if (previousProbeStructure != "")
	{
		previousProbeStructure = probeStructure->m_name;
	}
	else
	{
		previousProbeStructure = probeStructureName;
	}

	if (probeStructure != NULL)
	{
		delete(probeStructure);
	}


	probeStructure = new ProbeStructure(sceneName, probeStructureName);

	for (int i = 0; i < probeStructure->probeCount(); ++i)
	{
		SProbe gpuProbe;
		Probe* cpuProbe = probeStructure->getProbe(i);
		for (int j = 0; j < 9; j += 1)
		{
			gpuProbe.coefficients[(3 * j)] = cpuProbe->coeffs[j].x;
			gpuProbe.coefficients[(3 * j) + 1] = cpuProbe->coeffs[j].y;
			gpuProbe.coefficients[(3 * j) + 2] = cpuProbe->coeffs[j].z;

			gpuProbe.gradients[(9 * j) + 0] = cpuProbe->coeffGradients[j][0].x;
			gpuProbe.gradients[(9 * j) + 1] = cpuProbe->coeffGradients[j][0].y;
			gpuProbe.gradients[(9 * j) + 2] = cpuProbe->coeffGradients[j][0].z;

			gpuProbe.gradients[(9 * j) + 3] = cpuProbe->coeffGradients[j][1].x;
			gpuProbe.gradients[(9 * j) + 4] = cpuProbe->coeffGradients[j][1].y;
			gpuProbe.gradients[(9 * j) + 5] = cpuProbe->coeffGradients[j][1].z;

			gpuProbe.gradients[(9 * j) + 6] = cpuProbe->coeffGradients[j][2].x;
			gpuProbe.gradients[(9 * j) + 7] = cpuProbe->coeffGradients[j][2].y;
			gpuProbe.gradients[(9 * j) + 8] = cpuProbe->coeffGradients[j][2].z;

		}
		gpuProbe.position[0] = cpuProbe->position[0];
		gpuProbe.position[1] = cpuProbe->position[1];
		gpuProbe.position[2] = cpuProbe->position[2];
		gpuProbe.position[3] = 0.0f;
		probeList.probes[i] = gpuProbe;
	}

	for (int i = 0; i < probeStructure->m_dimensions.size(); ++i)
	{
		probeList.dimensions[i] = probeStructure->m_dimensions[i];
	}
	probeList.dimensions[3] = 0;;

	probeList.firstProbePosition[0] = probeStructure->m_firstProbePosition[0];
	probeList.firstProbePosition[1] = probeStructure->m_firstProbePosition[1];
	probeList.firstProbePosition[2] = probeStructure->m_firstProbePosition[2];
	probeList.firstProbePosition[3] = 0.0f;

	probeList.step = probeStructure->m_step;

	// Generate the SSBO holding probe information
	GLuint mySSBO = 0;
	glGenBuffers(1, &mySSBO);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, mySSBO);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(probeList), &probeList, GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, mySSBO);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	if (sampleSet)
	{
		sampleSet->probeStructure = probeStructure;
	}
}

void App::setSHTexturesUniforms(Args& args)
{
	args.setUniform("uSampler0", shTextures[0], Sampler());
	args.setUniform("uSampler1", shTextures[1], Sampler());
	args.setUniform("uSampler2", shTextures[2], Sampler());
	args.setUniform("uSampler3", shTextures[3], Sampler());
	args.setUniform("uSampler4", shTextures[4], Sampler());
	args.setUniform("uSampler5", shTextures[5], Sampler());
	args.setUniform("uSampler6", shTextures[6], Sampler());
	args.setUniform("uSampler7", shTextures[7], Sampler());
	args.setUniform("uSampler8", shTextures[8], Sampler());
}

void App::setProbeCoeffUniforms(Args& args, G3D::Array<G3D::Vector3>& coeffs)
{
	args.setUniform("coeff0", coeffs[0]);
	args.setUniform("coeff1", coeffs[1]);
	args.setUniform("coeff2", coeffs[2]);
	args.setUniform("coeff3", coeffs[3]);
	args.setUniform("coeff4", coeffs[4]);
	args.setUniform("coeff5", coeffs[5]);
	args.setUniform("coeff6", coeffs[6]);
	args.setUniform("coeff7", coeffs[7]);
	args.setUniform("coeff8", coeffs[8]);
}

void App::renderActors(RenderDevice* rd)
{
	Args args;

	G3D::Sampler shSampler = Sampler(G3D::WrapMode::CLAMP, G3D::InterpolateMode::NEAREST_NO_MIPMAP);
	
    //Iterate on each actor
	for (Actor& actor: actors)
	{
        //Don't draw if he's invisible
        if (!actor.isVisible) continue;
	
		// draw manipulator?
		actor.getManipulator()->render3D(rd);

        args = Args();
		args.setUniform("albedo", actor.albedo);
		args.setUniform("multiplier", shadingMultiplier);
		args.setUniform("r", 1.F);
		if (CPUInterpolation)
		{

			setProbeCoeffUniforms(args, actor.coefficients);

			//Accumulate the renders
			drawModel(rd, "SH_shader3.*", sphereModel, actor.getFrame(), args);    
		}
		else
		{
			
			args.setUniform("actorCentroidPosition", actor.getPosition());

			if (useSHGradients)
			{
				args.setMacro("USE_SH_GRADIENTS", 1);
			}

			args.setMacro("AMT_INTERP_PROBES", probeStructure ? probeStructure->probeCount() : 0);
			//args.setMacro("AMT_INTERP_PROBES", 9);
			//setSHTexturesUniforms(args);
			drawModel(rd, "SH_shader2.*", actor.getModel(), actor.getFrame(), args);
		}



    }

}

void App::drawProbes(RenderDevice* rd)
{
	G3D::Array<Probe*>& probeRenderList = showAllProbes ? probeStructure->probeList : probesToRender;
	Args args;

	for (Probe* probe : probeRenderList)
	{

		if (highlightProbes && !probe->bNeedsUpdate)
		{
			args.setUniform("multiplier", shadingMultiplier);
            args.setUniform("r", 1.f);

			setProbeCoeffUniforms(args, probe->coeffs);

			//Accumulate the renders
			drawModel(rd, "SH_shader3.*", sphereModel, probe->getPosition(), args);
		}
		else
		{
			if (useIrradianceTexture) args.setUniform("uSampler", Texture::fromFile("<white>"), Sampler());
			else args.setUniform("uSampler", probe->getTexture(0), Sampler());
			args.setUniform("multiplier", (float)1.0);
			drawModel(rd, "texture.*", sphereModel, probe->getPosition(), args);
			

		}


		if (bManipulateProbesEnabled)
		{
			probe->getManipulator()->render3D(rd);
		}
    }
}

void App::drawLights(RenderDevice* rd)
{
    Args args = Args();
    args.setUniform("multiplier", (float)1.0);

	for (SceneLight light: m_scene.m_sceneLights)
	{
		args.setUniform("uColor", light.color);
		drawModel(rd, "color.*", sphereModel, CFrame(light.position), args);
	}
}

void App::drawSurfaceSamples(RenderDevice* rd)
{

	int numOfsamplesToShow = atoi(maxSamplesPointsToDraw.c_str());

	for (int i =0; i < numOfsamplesToShow; i++){
		SceneSample& ss = sampleSet->m_samples[i];
		Color3 color = ss.irradiance;
		if (color == Color3(0, 0, 0))
		{
			if (!showDarkSamples) continue;
			color = Color3(1, 0, 0);
		}

		Draw::point(ss.position, rd, color * sampleMultiplier, 4.0f);

		if (showSampleNormals)
		{
			Draw::arrow(ss.position + (0.01 * ss.normal), ss.normal, rd, ss.irradiance, 0.1f); 
		}
	}

}

void App::drawScene(RenderDevice* rd)
{
    Args args;
    for (int i = 0; i < m_scene.m_models.size(); i++)
	{
        //if ((i == 5) && hideCeiling) continue;

        args = Args();
        args.setUniform("multiplier", sceneTextureIntensity);
		if (String(scenePane.selectedSceneList->selectedValue()) == "crytek-sponza")
		{
			//args.setUniform("uSampler", sceneTextures[0][i], Sampler());
			//drawModel(rd, "texture.*", sceneModels[i], CFrame(), args);
		}
		else
		{
			args.setUniform("uColor", m_scene.m_colors[i]);
			drawModel(rd, "color.*", m_scene.m_models[i], CFrame(), args);
		}
    }
}

void App::drawProbeLineSegments(RenderDevice* rd)
{
	for (int i = 0; i < probeStructure->lineArraySize(); i++){
        Draw::lineSegment(probeStructure->getLineSegment(i), rd, Color3::yellow());
     }

	for (int i = 0; i < probeStructure->probeCount(); i++){
		Probe* p = probeStructure->getProbe(i);
        Draw::lineSegment(LineSegment::fromTwoPoints(p->position, p->position+ p->normal*50.), rd, Color3::brown());
    }
}

/**
 * [App::onGraphics3D Test]
 * @param rd        [Le renderdevice]
 * @param surface3D [unused]
 */
void App::onGraphics3D(RenderDevice* rd, Array<shared_ptr<Surface> >& surface3D) {


	GBuffer::Specification gbufferSpec = m_gbufferSpecification;
	gbufferSpec.encoding[GBuffer::Field::WS_POSITION].format = ImageFormat::RGBA16F();
	extendGBufferSpecification(gbufferSpec);
	m_gbuffer->setSpecification(gbufferSpec);
	m_gbuffer->resize(m_framebuffer->width(), m_framebuffer->height());
	m_gbuffer->prepare(rd, activeCamera(), 0, -(float)previousSimTimeStep(), m_settings.hdrFramebuffer.depthGuardBandThickness, m_settings.hdrFramebuffer.colorGuardBandThickness);

	m_renderer->render(rd, m_framebuffer, scene()->lightingEnvironment().ambientOcclusionSettings.enabled ? m_depthPeelFramebuffer : shared_ptr<Framebuffer>(),
		scene()->lightingEnvironment(), m_gbuffer, surface3D);


	//if (submitToDisplayMode() == SubmitToDisplayMode::MAXIMIZE_THROUGHPUT) {
	//	swapBuffers();
	//}





 //   Args args;
 //   CFrame cframe;

    rd->setColorClearValue(Color3(0.3f, 0.3f, 0.3f));
    rd->setRenderMode(RenderDevice::RENDER_SOLID);
 //   
 //   //Set the framebuffer for drawing
    rd->pushState(m_framebuffer);

 //   //Clear the renderdevice
    //rd->clear();
	if (showSamples)
	{
		drawSurfaceSamples(rd);
	}

	//if (!hideActors)
	//{
	//	renderActors(rd);
	//}

	////Draw Probes
    if (probeStructure && (showInterpolationProbes || showAllProbes))
	{
		drawProbes(rd);
    }

	////Draw the lights
	//drawLights(rd);

 //   //Draw immobile scene
 //   if (!hideRoom)
	//{
	//	drawScene(rd);
	//}

	//if (CPUInterpolation)
	//{
	//	drawProbeLineSegments(rd);
	//	for (int i = 0; i < extrapolationTriangleVertices.size(); ++i)
	//	{
	//		Sphere s = Sphere(extrapolationTriangleVertices[i], 0.05);
	//		Draw::sphere(s, rd);
	//	}
	//}

 //   //Draw manipulator and other stuff
	//drawDebugShapes();
    rd->popState();

	//drawDebugShapes();
	rd->clear();
    m_film->exposeAndRender(rd, activeCamera()->filmSettings(), m_framebuffer->texture(0), 0, 0);
 //   timer++;

	//// I'm only doing this because I can't get rd->screenshot() to save to a good filename
	//if (numPassesLeft > 0)
	//{
	//	const shared_ptr<Image> screen(rd->screenshotPic());

	//	FileSystem::ListSettings ls;
	//	ls.includeParentPath = false;
	//	ls.recursive = false;
	//	G3D::Array<G3D::String> screenshotList;
	//	FileSystem::list("C:/temp/CurrentOptimization/screens/*", screenshotList, ls);
	//	int num = screenshotList.size();

	//	String filename = "C:/temp/CurrentOptimization/screens/" + String(std::to_string(num).c_str()) + ".jpg";
	//	debugPrintf("%s", filename.c_str());
	//	screen->save(filename);

	//	debugPrintf("%d optimization passes left \n", numPassesLeft);
	//	numPassesLeft--;
	//}
}


void App::drawModel(RenderDevice* rd, String shaderName, shared_ptr<ArticulatedModel> model, CFrame frame, Args args){
    args.setUniform("intensity", 0.1f); 

    Array<shared_ptr<Surface> > surfaceArray = Array<shared_ptr<Surface> >();
    model->pose(surfaceArray, frame);
	for (int i = 0; i < surfaceArray.size(); i++){
		const shared_ptr<UniversalSurface> &surface = dynamic_pointer_cast<UniversalSurface>(surfaceArray[i]);
		CFrame cframe;
		if (notNull(surface)) {

			/*const shared_ptr<Texture> tex = surface->material()->bsdf()->lambertian().texture();
			args.setUniform("uSampler", tex, Sampler());*/

			surface->getCoordinateFrame(cframe);
			args.setUniform("MVP", rd->projectionMatrix() * (rd->cameraToWorldMatrix().inverse() * frame));
			rd->setObjectToWorldMatrix(cframe);
			surface->gpuGeom()->setShaderArgs(args);
			if (shaderName.compare("color.*") == 0){
				LAUNCH_SHADER("../data-files/Shaders/color.*", args);    
			} 
			else if (shaderName.compare("texture.*") == 0){
				LAUNCH_SHADER("../data-files/Shaders/texture.*", args);    
			} 
			else if (shaderName.compare("depth.*") == 0){
				LAUNCH_SHADER("../data-files/Shaders/depth.*", args);    
			} 
			else if (shaderName.compare("SH_shader2.*") == 0){
				LAUNCH_SHADER("../data-files/Shaders/SHshader2.*", args);    
			} 
			else if (shaderName.compare("SH_shader3.*") == 0){
				LAUNCH_SHADER("../data-files/Shaders/SH_shader3.*", args);    
			} 
			else if (shaderName.compare("SHFunctions.*") == 0){
				LAUNCH_SHADER("../data-files/Shaders/SHFunctions.*", args);    
			} 
			else{
				LAUNCH_SHADER("../data-files/Shaders/SH_shader.*", args);    
			}
		}
	}
}

void App::addActor(String name, String filename, Point3 position, float scale, shared_ptr<Texture> texture){
    ArticulatedModel::Specification spec;
    spec.filename       = System::findDataFile(filename);
    spec.scale          = scale;
    Actor myActor = Actor(name, spec, position, texture);

    actors.append(myActor);
    //addWidget(myActor.getManipulator());

}

void App::addActor(String name, shared_ptr<ArticulatedModel> model, Point3 position, float scale, shared_ptr<Texture> texture, bool useManipulator, Vector3 albedo)
{
	Actor myActor = Actor(name, model,position, texture, true);
	myActor.albedo = albedo;
    actors.append(myActor);

	if (useManipulator)
	{
	    addWidget(myActor.getManipulator());
	}
	
}



void App::addModel(String filename, Color3 color)
{
	//if (!//loadOBJ(filename))
	//{
	//	debugPrintf("Couldn't load file %s", filename);
	//	return;
	//}
	//m_scene.m_colors.push_back(Vector3(color));
}

void App::addLight(Point3 pos, Color3 color)
{
	SceneLight light;
	light.position = pos;
	light.color = color;
	sceneLights.append(light);
}

void App::addLight(SceneLight light)
{
	sceneLights.append(light);
}

void App::debugStringArray(Array<String> myArray){
    for (int i =0; i < myArray.size(); i++){
        debugPrintf("array[%d] = %s ",i,myArray[i].c_str());
    }
    debugPrintf("\n");
}

void App::debugIntArray(Array<int> myArray){
    for (int i =0; i < myArray.size(); i++){
        debugPrintf("array[%d] = %d ",i,myArray[i]);
    }
    debugPrintf("\n");
}
void App::debugFloatArray(Array<float> myArray){
    for (int i =0; i < myArray.size(); i++){
        debugPrintf("array[%d] = %f ",i,myArray[i]);
    }
    debugPrintf("\n");
}

void App::loadSHTextures(){
    debugPrintf("Loading SH Textures...\n");

    for (int i = 0; i < 100; i++){
        char str[8];
        sprintf(str, "%d", i);
        String filename = "../data-files/SH/512x256/SHRender";
        filename.append(str);
        filename.append(".png");
        shTextures.append(Texture::fromFile(filename));
    }
}

void App::makeGui() {
    shared_ptr<GuiWindow> gui = GuiWindow::create("Parameters");    
    GuiPane* pane = gui->pane();
	    
	GuiTabPane* tabPane = pane->addTabPane();
	GuiPane* tab = tabPane->addTab("General Controls");

    tab->beginRow();
	tab->addCheckBox("Interpolate Coefficients", &interpolateCoefficients);
	tab->addButton("displace", GuiControl::Callback(this, &App::displaceProbes), GuiTheme::TOOL_BUTTON_STYLE);

    tab->endRow();

    tab->beginRow();
	tab->addButton("Sample Scene", GuiControl::Callback(this, &App::computeSceneSamples),GuiTheme::TOOL_BUTTON_STYLE);
	tab->addButton("Compute samplesRGB", GuiControl::Callback(this, &App::computeSamplesRGB), GuiTheme::TOOL_BUTTON_STYLE);
	tab->addButton("Compute samplesRGB_ref", GuiControl::Callback(this, &App::computeSamplesRGBRef), GuiTheme::TOOL_BUTTON_STYLE);
	tab->addButton("Compute triplets", GuiControl::Callback(this, &App::computeTriplets), GuiTheme::TOOL_BUTTON_STYLE);
	tab->addTextBox("samplesToSave", samplesToSave);
    tab->addCheckBox("write samples", &saveSample);
    tab->addCheckBox("MATLAB", &useMatlabOptimization);
    tab->addButton("Optimization pass", GuiControl::Callback(this, &App::startOptimizationPasses),GuiTheme::TOOL_BUTTON_STYLE);
	tab->addTextBox("Num", &tbNumPassesLeft);
    tab->addCheckBox("log", &logSampleSet);
    tab->endRow();

    tab->beginRow();
	tab->addTextBox("NumSamples", &maxSamplesPointsToDraw);
	tab->addCheckBox("Show", &showSamples);
	tab->addCheckBox("Show normals", &showSampleNormals);
	tab->addCheckBox("Show dark samples", &showDarkSamples);
	tab->addSlider("F", &sampleMultiplier, 1.0f, 5.f);
    tab->addButton("Find Initial Conditions", GuiControl::Callback(this, &App::findBestInitialConditions),GuiTheme::TOOL_BUTTON_STYLE);
	tab->addTextBox("NumTries", &m_sNumICTries);
	tab->addTextBox("NumProbes", &m_sNumICProbes);
	tab->endRow();

	tab = tabPane->addTab("Rendering");
	tab->beginRow();
	tab->addCheckBox("Direct", &(bRenderDirect));
	tab->addCheckBox("Indirect (probes)", &(bRenderIndirect));
	tab->addCheckBox("* BRDF", &(bRenderMultiplyIndirectByBRDF));
	tab->addSlider("Mutiplier", &shadingMultiplier, 0.0f, 5.0f);
	tab->addCheckBox("AOF", &(bRenderAO));
	tab->addCheckBox("Shadow Maps", &(bRenderShadowMaps));
	tab->endRow();

	tab->beginRow();
	tab->addCheckBox("Show Interp Probes", &showInterpolationProbes);
	tab->addCheckBox("Show ALL Probes", &showAllProbes);
	tab->addCheckBox("CPUInterpolation", &CPUInterpolation);
	tab->addCheckBox("highlight probes", &highlightProbes);
	tab->addButton("OfflineRender", GuiControl::Callback(this, &App::offlineRender), GuiTheme::TOOL_BUTTON_STYLE);
	tab->endRow();

	tab = tabPane->addTab("Actors");
	tab->beginRow();
	tab->addButton("Add sphere", GuiControl::Callback(this, &App::addOneActor), GuiTheme::TOOL_BUTTON_STYLE);
	tab->addButton("Add sceneActor", GuiControl::Callback(this, &App::addOneActorSq), GuiTheme::TOOL_BUTTON_STYLE);
	tab->addButton("Clear actors", GuiControl::Callback(this, &App::clearAllActors), GuiTheme::TOOL_BUTTON_STYLE);
	tab->addTextBox("X", &actorSpawnX);
	tab->addTextBox("Y", &actorSpawnY);
	tab->addTextBox("Z", &actorSpawnZ);
	tab->endRow();
    for (int i = 0 ; i < actors.size(); i++)
	{
		tab->beginRow();
		tab->addLabel(GuiText(actors[i].name));
		tab->addCheckBox("visible", &actors[i].updateVisibility);
		tab->addSlider("Phong",    &actors[i].phongExponent, 1.0f, 60.f);
		tab->addSlider("Max band",    &actors[i].drawBand, 0, maxDrawBand);
		tab->addCheckBox("show this band only", &actors[i].showOnlyOneBand);
		tab->endRow();
    }

	addScenePane(tabPane);
    gui->pack();
    addWidget(gui);
    gui->moveTo(Point2(10, 50));
}

void App::addScenePane(GuiTabPane* tabPane)
{
	GuiPane* tab = tabPane->addTab("Scene");

	G3D::Array<G3D::String> sceneList;

	FileSystem::ListSettings ls;
	ls.includeParentPath = false;
	ls.recursive = false;
	
	FileSystem::list("../data-files/Scenes/*", sceneList, ls);

	scenePane.selectedSceneList = tab->addDropDownList("Scene", sceneList, NULL, GuiControl::Callback(this, &App::updateSelectedScenePane));
	
	G3D::String selectedScene = scenePane.selectedSceneList->selectedValue();

	tab->beginRow();
	G3D::Array<G3D::String> probeStructureList;
	FileSystem::list("../data-files/Scenes/" + selectedScene + "/ProbeStructures/*", probeStructureList, ls);
	scenePane.probeStructureList = tab->addDropDownList("ProbeStructure", probeStructureList, NULL, GuiControl::Callback(this, &App::updateProbeStructure));
	tab->addButton(GuiText("Switch Back"), GuiControl::Callback(this, &App::loadPreviousProbeStructure), GuiTheme::TOOL_BUTTON_STYLE );
	tab->addButton(GuiText("Edit mode"), GuiControl::Callback(this, &App::switchEditProbeStructure), GuiTheme::TOOL_BUTTON_STYLE);
	tab->addButton(GuiText("Save"), GuiControl::Callback(this, &App::saveProbeStructureUpdate), GuiTheme::TOOL_BUTTON_STYLE);
	tab->addButton(GuiText("Update All"), GuiControl::Callback(this, &App::saveProbeStructureUpdateAll), GuiTheme::TOOL_BUTTON_STYLE );
	
	tab->endRow();


	G3D::Array<G3D::String> sampleSetList;
	FileSystem::list("../data-files/Scenes/" + selectedScene + "/SampleSets/*", sampleSetList, ls);
	scenePane.sampleSetList = tab->addDropDownList("SampleSet", sampleSetList, NULL, GuiControl::Callback(this, &App::updateSampleSet));
	tab->addButton(GuiText("Clear"), GuiControl::Callback(this->sampleSet, &SceneSampleSet::clear), GuiTheme::TOOL_BUTTON_STYLE);
}

void App::saveProbeStructureUpdate()
{
	probeStructure->updateProbes(false);
	updateProbeStructure();
}

void App::saveProbeStructureUpdateAll()
{
	probeStructure->updateProbes(true);
	updateProbeStructure();
}

void App::switchEditProbeStructure()
{
	if (!bManipulateProbesEnabled)
	{
		for (Probe* p : probeStructure->probeList)
		{
			addWidget(p->getManipulator());
		}
		bManipulateProbesEnabled = true;
	}
	else
	{
		for (Probe* p : probeStructure->probeList)
		{
			removeWidget(p->getManipulator());
		}
		bManipulateProbesEnabled = false;
	}
	
}

void App::updateSelectedScenePane()
{
	G3D::String selectedScene = scenePane.selectedSceneList->selectedValue();
	debugPrintf("update ldofgio");

	loadScene(selectedScene);

	FileSystem::ListSettings ls;
	ls.includeParentPath = false;
	ls.recursive = false;

	G3D::Array<G3D::String> probeStructureList;
	FileSystem::list("../data-files/Scenes/" + selectedScene + "/ProbeStructures/*", probeStructureList, ls);
	scenePane.probeStructureList->setList(probeStructureList);

	G3D::Array<G3D::String> sampleSetList;
	FileSystem::list("../data-files/Scenes/" + selectedScene + "/SampleSets/*", sampleSetList, ls);
	scenePane.sampleSetList->setList(sampleSetList);
}


void App::updateProbeStructure()
{
	G3D::String selectedScene = scenePane.selectedSceneList->selectedValue();
	G3D::String probeStructureName = scenePane.probeStructureList->selectedValue();
	initializeProbeStructure(selectedScene, probeStructureName);
}

void App::loadPreviousProbeStructure()
{
	if (previousProbeStructure != "")
	{
		initializeProbeStructure(scenePane.selectedSceneList->selectedValue(), previousProbeStructure);
	}
	if (sampleSet)
	{
		sampleSet->probeStructure = probeStructure;
	}
}

void App::updateSampleSet()
{
	G3D::String selectedScene = scenePane.selectedSceneList->selectedValue();
	G3D::String sampleSetName = scenePane.sampleSetList->selectedValue();
	
	//loadSurfaceSamples(selectedScene, sampleSetName);
	sampleSet = new SceneSampleSet(std::string(selectedScene.c_str()), std::string(sampleSetName.c_str()), m_scene.m_scale);
	if (probeStructure)
	{
		sampleSet->probeStructure = probeStructure;
	}
}

void App::updateSelectedSceneTextures(){
	debugPrintf("Wow! The dropdown list index has been set to %d!!!!\n", sampleDropDownIndex);
}


void App::saveCoeffs(String path){
    std::fstream coeffFile;
	debugPrintf("Saving coeffs \n");

    coeffFile.open(path.c_str(), std::fstream::out);
	if (!coeffFile)
	{
		throw std::runtime_error("..");
	}
    for (int i = 0; i < actors[0].coefficients.size(); i++){
			debugPrintf("%f\n",actors[0].coefficients[i].x);
            coeffFile<<actors[0].coefficients[i].x << "\n";
            coeffFile<<actors[0].coefficients[i].y << "\n";
            coeffFile<<actors[0].coefficients[i].z << "\n";
        
    }
    coeffFile.close();
}

void App::loadSampleCoeffs()
{
	std::fstream coeffFile;
	debugPrintf("loading sampled coeffs\n");
    coeffFile.open("../data-files/ProbeStructures/computed/GeneratedCoeffFromSamples.txt", std::fstream::in);

	Array<Vector3> loadedCoeffs;

	for (int i = 0; i < 9; i++)
	{
		std::string line;
		std::getline(coeffFile, line);
		float cr = std::stof(line.c_str());

		std::getline(coeffFile, line);
		float cg = std::stof(line.c_str());

		std::getline(coeffFile, line);
		float cb = std::stof(line.c_str());

		Vector3 v = Vector3(cr, cg, cb);

		loadedCoeffs.append(v);
	}
	debugPrintf("%s\n",printVector3Array(loadedCoeffs).c_str());
	actors[0].coefficients = loadedCoeffs;
}

void App::onAI(){
    GApp::onAI();

	if (actors.size() > 0)
	{
		screenPrintf("Actor pos: %s", actors[0].getPosition().toString().c_str());

	}

	if (numPassesLeft > 0)
	{
		tryOptimization();
	}

	if (probeStructure)
	{
		screenPrintf("ProbeStructure: %s", probeStructure->m_name.c_str());
	}
	
	//TODO: probe ray casting one day? if i find a real use for it

    for (int i = 0 ; i < actors.size(); i++){


		if (interpolateCoefficients)
		{
			bool keepProbes = (i == 0);
			//actors[i].coefficients = triLinearInterpolation(i, &probeIndices, &weights, keepProbes);

			ProbeInterpolationRecord iRec = probeStructure->getInterpolationProbeIndicesAndWeights(actors[i].getPosition());

			//if (tetrahedralInterpolation(actors[i], &probeIndices, &weights))
			{
				actors[i].coefficients = getInterpolatedCoeffs(iRec, actors[i].drawBand);
			}
			
			if (showInterpolationProbes)
			{
				probesToRender.clear();
				for (int i = 0; i < iRec.probeIndices.size(); ++i)
				{
					probesToRender.push_back(probeStructure->getProbe(iRec.probeIndices[i]));
				}
			}


		}
        if (actors[i].updateVisibility)
		{
            if (actors[i].isVisible)
			{
                actors[i].isVisible = false;
                removeWidget(actors[i].getManipulator());
            }
            else
			{
                actors[i].isVisible = true;
                addWidget(actors[i].getManipulator());
            }
            actors[i].updateVisibility = false;
        }
    }


}

String App::printVector3Array(Array<Vector3> arr){
    String s = "";
    char str[50];
    for (int i = 0; i < arr.size(); i++){
        sprintf(str, "[%d] x: %f, y: %f, z: %f ",i, arr[i].x, arr[i].y, arr[i].z);
        s.append(str);
    }

    return s;
}

void App::loadCornell()
{
    ArticulatedModel::Specification spec = ArticulatedModel::Specification();
    m_scene.m_meshShapes = Array<MeshShape*>();
	m_scene.m_scale = 0.01f;
    //loadOBJ("../data-files/Scenes/zcbox/meshes/cbox_back.obj",		m_scene.scale());
    //loadOBJ("../data-files/Scenes/zcbox/meshes/cbox_ceiling.obj",	m_scene.scale());
    //loadOBJ("../data-files/Scenes/zcbox/meshes/cbox_floor.obj",		m_scene.scale());
    //loadOBJ("../data-files/Scenes/zcbox/meshes/cbox_greenwall.obj", m_scene.scale());
    //loadOBJ("../data-files/Scenes/zcbox/meshes/cbox_largebox.obj",	m_scene.scale());
    //loadOBJ("../data-files/Scenes/zcbox/meshes/cbox_luminaire.obj", m_scene.scale());
    //loadOBJ("../data-files/Scenes/zcbox/meshes/cbox_redwall.obj",	m_scene.scale());
    //loadOBJ("../data-files/Scenes/zcbox/meshes/cbox_smallbox.obj",	m_scene.scale());

    loadCornellColors();
}

void App::loadCornellColors()
{
    m_scene.m_colors.append(Vector3(1,1,1));
	m_scene.m_colors.append(Vector3(1,1,1));
	m_scene.m_colors.append(Vector3(1,1,1));
	m_scene.m_colors.append(Vector3(0,1,0));
	m_scene.m_colors.append(Vector3(0,0,1));
	m_scene.m_colors.append(Vector3(1,1,1));
	m_scene.m_colors.append(Vector3(1,0,0));
	m_scene.m_colors.append(Vector3(0,0,1));

}


void App::loadSponza(){
	m_scene.m_scale = 1.0f;
    //loadOBJ("../data-files/Scenes/crytek-sponza/sponza_v2_scaled.obj");

    makeSponzaTextures();
	addLight(Point3(0,3,0), Color3(1,0,0));
}

//comment my friend
void App::loadZScene(){
	m_scene.m_scale = 1.0f;
    //loadOBJ("../data-files/Z2/objparts/white_connectedtored.obj");
    //loadOBJ("../data-files/Z2/objparts/bluewall_connectedtogreen.obj");
    //loadOBJ("../data-files/Z2/objparts/white_connectedtogreen.obj");
    //loadOBJ("../data-files/Z2/objparts/bluewall_connectedtored.obj");
    //loadOBJ("../data-files/Z2/objparts/floor.obj");
    //loadOBJ("../data-files/Z2/objparts/ceiling.obj");
    //loadOBJ("../data-files/Z2/objparts/largeredwall.obj");
    //loadOBJ("../data-files/Z2/objparts/smallredwall.obj");
    //loadOBJ("../data-files/Z2/objparts/largegreenwall.obj");
    //loadOBJ("../data-files/Z2/objparts/smallgreenwall2.obj");

    makeZTextures();

	addLight(Point3(0,4,0), Color3(1,1,1));
	addLight(Point3(-3.5,4,-3.5), Color3(1,1,1));
	addLight(Point3(3.5,4,3.5), Color3(1,1,1));
}

void App::makeSponzaTextures(){
	sceneTextures = Array<Array<shared_ptr<Texture> > >();
	sceneTextures.append(Array<shared_ptr<Texture> >());

	sceneTextures[0].append(Texture::fromFile("../data-files/Scenes/crytek-sponza/sponz_diff2.tga"));

}

// LLoad the scene textures
void App::makeZTextures(){
	sceneTextures = Array<Array<shared_ptr<Texture> > >();
	sceneTextures.append(Array<shared_ptr<Texture> >());

	Array<String> textures;
	FileSystem::list("../data-files/Scenes/ZScene/cycles_renders/*", textures);
	String s = System::findDataFile("../data-files/Scenes/ZScene/cycles_renders/white_redconnect_3000.png");

    //sceneTextures[0].append(Texture::fromFile(textures[0]));
    //sceneTextures[0].append(Texture::fromFile("/Scenes/ZScene/cycles_renders/blue_greenside_3000_rotated.png"));
    //sceneTextures[0].append(Texture::fromFile("/Scenes/ZScene/cycles_renders/white_greenconnect_3000.png"));
    //sceneTextures[0].append(Texture::fromFile("/Scenes/ZScene/cycles_renders/blue_redside_3000_rotated.png"));
    //sceneTextures[0].append(Texture::fromFile("/Scenes/ZScene/cycles_renders/floor_3000.png"));
    //sceneTextures[0].append(Texture::fromFile("/Scenes/ZScene/cycles_renders/ceiling_3000.png"));
    //sceneTextures[0].append(Texture::fromFile("/Scenes/ZScene/cycles_renders/big_red_3000.png"));
    //sceneTextures[0].append(Texture::fromFile("/Scenes/ZScene/cycles_renders/small_red_3000_rotated.png"));
    //sceneTextures[0].append(Texture::fromFile("/Scenes/ZScene/cycles_renders/big_green_3000.png"));
    //sceneTextures[0].append(Texture::fromFile("/Scenes/ZScene/cycles_renders/small_green_3000.png"));

	m_scene.m_colors.append(Vector3(1,1,1));
	m_scene.m_colors.append(Vector3(0,0,1));
	m_scene.m_colors.append(Vector3(1,1,1));
	m_scene.m_colors.append(Vector3(0,0,1));
	m_scene.m_colors.append(Vector3(1,1,1));
	m_scene.m_colors.append(Vector3(1,1,1));
	m_scene.m_colors.append(Vector3(1,0,0));
	m_scene.m_colors.append(Vector3(1,0,0));
	m_scene.m_colors.append(Vector3(0,1,0));
	m_scene.m_colors.append(Vector3(0,1,0));

	//System::findDataFile("../data-files/objs/sphere.obj");
}



void App::offlineRender()
{
	SOfflineRenderingOptions offlineRenderingOptions;
	offlineRenderingOptions.numSamples = "256";
	offlineRenderingOptions.height = "256";
	offlineRenderingOptions.width = "512";
	std::stringstream command;
	command << "cmd /c \"cd C:\\git\\AutoProbePlacement\\AutoProbePlacement\\data-files\\scripts";
	command << " && C:\\Users\\Joel\\Anaconda2\\python.exe CreateSceneRender.py zcbox 1Probe Probes ";

	command << m_debugCamera->frame().translation.x << " ";
	command << m_debugCamera->frame().translation.y << " ";
	command << m_debugCamera->frame().translation.z << " ";
	command << m_debugCamera->frame().lookVector().x << " ";
	command << m_debugCamera->frame().lookVector().y << " ";
	command << m_debugCamera->frame().lookVector().z << " ";
	command << actors[0].getPosition().x << " ";
	command << actors[0].getPosition().y << " ";
	command << actors[0].getPosition().z << " ";
	command << offlineRenderingOptions.numSamples.c_str() << " ";
	command << offlineRenderingOptions.width.c_str() << " ";
	command << offlineRenderingOptions.height.c_str() << " ";
	command << "\"";

	debugPrintf("%s\n", command.str().c_str());
	runCommand(command.str());
}