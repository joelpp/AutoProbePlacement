/** \file App.cpp */
#include <fstream>

#include "App.h"
#include "SceneSampleSet.h"
#include "ProbeRenderer.h"
#include "Helpers.h"
#include "Integrator.h"

/*
    Debugging Helpers
*/

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
    
	createNewSampleSetWindow();
	createNewProbeWindow();
    createNewProbeStructureWindow();
    createCopyProbeStructureWindow();
    createRenameOptimizationWindow();

	probeStructurePane = NULL;	
	m_probeStructure = NULL;
    m_scene = NULL;
	sampleSet = NULL;

	integratorList.resize(EIntegrator::NUM_INTEGRATORS);
	integratorList[EIntegrator::Path] = "path";
	integratorList[EIntegrator::Indirect] = "path_samples";
	integratorList[EIntegrator::Direct] = "direct";

	filmTypeList.resize(EFilmType::NUM_FILM_TYPES);
	filmTypeList[EFilmType::HDR] = "hdrfilm";
	filmTypeList[EFilmType::LDR] = "ldrfilm";

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
    spec.scale = 1.0f;
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
	saveSample = true;
	useBakedSceneTextures = true;
	interpolateCoefficients = false;
	logSampleSet = false;
	showSampleNormals = false;
	showDarkSamples = true;
	showSamples = false;
	useSHGradients = false;
	bManipulateProbesEnabled = false;
	hideActors = false;
    bOptimizeForCoeffs = false;
	bRenderDirect = true;
	bRenderIndirect = true;
	bRenderMultiplyIndirectByBRDF = false;
    bShouldUpdateProbeStructurePane = true;
    bKeepRefValuesOnNewOptimization = false;
    bTakeRefScreenshot = false;
	bOptimizeWithMitsubaSamples = false;
    //Decide how many bands we want to use for the interpolation
    numBands = 2;
    totalNumberOfCoeffs = 0;
    maxDrawBand = 8;
    shadingMultiplier = 1.0f;
    sampleMultiplier = 4.0f;
    lightIntensity = 1;
    sceneTextureIntensity = 1.0;
    extrapolationT = 0;
    smallestBaryCoord = -1;
    numPassesLeft = 0;
    timer = 0;
	maxProbeStepLength = "0.1";

    numOptimizationSamples = "0";
    maxSamplesPointsToDraw = "0";
    samplesToSave = new String("0");
    SampleSetOutputName = new String("SampleSetOutputName");
    actorSpawnX = String("0");
    actorSpawnY = String("0");
    actorSpawnZ = String("0");
    tbNumPassesLeft = String("1");
	shadingSHBand = "9";
	optimizationSHBand = "9";
	m_activeCamera->setPosition(Point3(0,3,3));
	sampleDropDownIndex = 0;
	tetgenString = "../data-files/tetgen/44/";
    previousProbeStructure = "";

    for (int i = 0; i <= numBands; i++)
	{
		totalNumberOfCoeffs += 2*i + 1; 
	}

    createDeveloperHUD();

    //Add sliders for rendering parameters
	addOneActor();

    //Load the textures containing the SH coefficient values
    //loadSHTextures();

	offlineRenderingOptions.numSamples = "256";
	offlineRenderingOptions.height = "256";
	offlineRenderingOptions.width = "512";
	offlineRenderingOptions.gamma = "2.2";

	setActiveCamera(m_debugCamera);

	loadOptions();
	makeGui();
}//end of onInit 

void App::loadScene(String sceneName)
{
	lightPositions.clear();

	m_scene = new JScene(sceneName);

    String basePath = "C:\\git\\AutoProbePlacement\\AutoProbePlacement\\data-files\\Scenes\\" + sceneName + "\\";
	try
	{
		GApp::loadScene(basePath + sceneName + ".Scene.Any");
	}
	catch (std::exception e)
	{
		throw e;
	}
    setActiveCamera(m_debugCamera);

    m_currentOptimization = folderCount(basePath + "Optimizations\\") - 1;
    debugPrintf("%d\n", m_currentOptimization);

}

void App::initializeProbeStructure(String sceneName, String probeStructureName)
{
	if (previousProbeStructure != "")
	{
		previousProbeStructure = m_probeStructure->m_name;
	}
	else
	{
		previousProbeStructure = probeStructureName;
	}

	if (m_probeStructure != NULL)
	{
		delete(m_probeStructure);
	}

	try
	{
		m_probeStructure = new ProbeStructure(sceneName, probeStructureName);
	}
	catch (std::exception e)
	{
		throw e;
	}

	if (sampleSet)
	{
		sampleSet->probeStructure = m_probeStructure;
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
		//actor.getManipulator()->render3D(rd);

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

			args.setMacro("AMT_INTERP_PROBES", m_probeStructure != NULL ? m_probeStructure->probeCount() : 0);
			//args.setMacro("AMT_INTERP_PROBES", 9);
			//setSHTexturesUniforms(args);
			drawModel(rd, "SH_shader2.*", actor.getModel(), actor.getFrame(), args);
		}
    }
}

void App::drawProbes(RenderDevice* rd)
{
	G3D::Array<Probe*>& probeRenderList = showAllProbes ? m_probeStructure->probeList : probesToRender;

	Args args;
	args.setUniform("multiplier", shadingMultiplier);
	args.setUniform("r", 1.f);

	if (bManipulateProbesEnabled && m_probeStructure->eType() == EProbeStructureType::Trilinear)
	{
		Probe* p = m_probeStructure->getProbe(0);
		setProbeCoeffUniforms(args, p->coeffs);

		drawModel(rd, "SH_shader3.*", sphereModel, p->getPosition(), args);

		//float step = m_probeStructure->m_step;
		float step = std::stof(probeStructurePanelOptions.step.c_str());

		Color3 white = Color3::white();
		shared_ptr<Texture> whiteTexture = Texture::white();
		Vector3 manipulatorPos = p->getManipulator()->frame().translation;
		Array<Vector3> positions;
		Sampler sampler = Sampler();
		for (int i = 0; i < m_probeStructure->m_dimensions[0]; ++i)
		{
			for (int j = 0; j < m_probeStructure->m_dimensions[1]; ++j)
			{
				for (int k = 0; k < m_probeStructure->m_dimensions[2]; ++k)
				{
					Vector3 centerPos = manipulatorPos + Vector3(i * step, j * step, k * step);
					//positions.push_back(manipulatorPos + Vector3(i * step, j * step, k * step)); 
					//Draw::sphere(Sphere(centerPos, 0.1f), rd, white);
					args.setUniform("uSampler", whiteTexture, sampler);
					drawModel(rd, "texture.*", sphereModel, centerPos, args);
				}
			}
		}
		//Draw::points(positions, rd, white, sampleMultiplier);

	}

	for (Probe* probe : probeRenderList)
	{

		if (highlightProbes && !probe->bNeedsUpdate)
		{


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

		//if (bManipulateProbesEnabled)
		//{
		//	probe->getManipulator()->render3D(rd);
		//}
    }
}

void App::drawLights(RenderDevice* rd)
{
    Args args = Args();
    args.setUniform("multiplier", (float)1.0);

	for (SceneLight light: m_scene->m_sceneLights)
	{
		args.setUniform("uColor", light.color);
		drawModel(rd, "color.*", sphereModel, CFrame(light.position), args);
	}
}

void App::drawSurfaceSamples(RenderDevice* rd)
{
	int numOfsamplesToShow = atoi(maxSamplesPointsToDraw.c_str());

	Draw::points(sampleSet->m_points, rd, sampleSet->m_colors, sampleMultiplier);
}

void App::drawScene(RenderDevice* rd)
{
    Args args;
    for (int i = 0; i < m_scene->m_models.size(); i++)
	{
        //if ((i == 5) && hideCeiling) continue;

        args = Args();
        args.setUniform("multiplier", sceneTextureIntensity);
		if (selectedSceneName() == "crytek-sponza")
		{
			//args.setUniform("uSampler", sceneTextures[0][i], Sampler());
			//drawModel(rd, "texture.*", sceneModels[i], CFrame(), args);
		}
		else
		{
			args.setUniform("uColor", m_scene->m_colors[i]);
			drawModel(rd, "color.*", m_scene->m_models[i], CFrame(), args);
		}
    }
}

void App::drawProbeLineSegments(RenderDevice* rd)
{
	for (int i = 0; i < m_probeStructure->lineArraySize(); i++){
        Draw::lineSegment(m_probeStructure->getLineSegment(i), rd, Color3::yellow());
     }

	for (int i = 0; i < m_probeStructure->probeCount(); i++){
		Probe* p = m_probeStructure->getProbe(i);
        Draw::lineSegment(LineSegment::fromTwoPoints(p->position, p->position+ p->normal*50.), rd, Color3::brown());
    }
}


/**
 * [App::onGraphics3D Test]
 * @param rd        [Le renderdevice]
 * @param surface3D [unused]
 */
void App::onGraphics3D(RenderDevice* rd, Array<shared_ptr<Surface> >& surface3D) {
	rd->setColorClearValue(Color3(0.3f, 0.3f, 0.3f));
	//rd->clear();
	GBuffer::Specification gbufferSpec = m_gbufferSpecification;
	gbufferSpec.encoding[GBuffer::Field::WS_POSITION].format = ImageFormat::RGBA16F();
	extendGBufferSpecification(gbufferSpec);
	m_gbuffer->setSpecification(gbufferSpec);
	m_gbuffer->resize(m_framebuffer->width(), m_framebuffer->height());
	m_gbuffer->prepare(rd, activeCamera(), 0, -(float)previousSimTimeStep(), m_settings.hdrFramebuffer.depthGuardBandThickness, m_settings.hdrFramebuffer.colorGuardBandThickness);

	m_renderer->render(rd, m_framebuffer, scene()->lightingEnvironment().ambientOcclusionSettings.enabled ? m_depthPeelFramebuffer : shared_ptr<Framebuffer>(),
		scene()->lightingEnvironment(), m_gbuffer, surface3D);

    //rd->setRenderMode(RenderDevice::RENDER_SOLID);
 //   
 //   //Set the framebuffer for drawing
    rd->pushState(m_framebuffer);

	if (showSamples)
	{
		drawSurfaceSamples(rd);
	}

	////Draw Probes
    if (m_probeStructure && (showInterpolationProbes || showAllProbes))
	{
		drawProbes(rd);
    }

    renderActors(rd);

 //   //Draw manipulator and other stuff
    rd->popState();

	//rd->clear();
    m_film->exposeAndRender(rd, activeCamera()->filmSettings(), m_framebuffer->texture(0), 0, 0);

    if (optimizing || numPassesLeft > 0)
    {
    	const shared_ptr<Image> screen(rd->screenshotPic());

    	FileSystem::ListSettings ls;
    	ls.includeParentPath = false;
    	ls.recursive = false;
    	G3D::Array<G3D::String> screenshotList;
    	FileSystem::list(currentOptimizationFolderPath() + "/screens/*", screenshotList, ls);
    	int num = screenshotList.size();

    	String filename = currentOptimizationFolderPath() + "/screens/" + String(std::to_string(num).c_str()) + ".jpg";
    	debugPrintf("%s", filename.c_str());
    	screen->save(filename);

    	debugPrintf("%d optimization passes left \n", numPassesLeft);
    }

    if (bTakeRefScreenshot)
    {
        String filename = currentOptimizationFolderPath() + "/ref_probes.jpg";

        const shared_ptr<Image> screen(rd->screenshotPic());
        screen->save(filename);

        bTakeRefScreenshot = false;
    }
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
	    //addWidget(myActor.getManipulator());
	}
	
}

void App::addModel(String filename, Color3 color)
{
	//if (!//loadOBJ(filename))
	//{
	//	debugPrintf("Couldn't load file %s", filename);
	//	return;
	//}
	//m_scene->m_colors.push_back(Vector3(color));
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

void App::clearSampleSetValues()
{
	if (sampleSet != NULL)
	{
sampleSet->clearValues();
    }
}

void App::clearSampleSetPositions()
{
    if (sampleSet != NULL)
    {
        sampleSet->clearPositions();
    }
}

void App::reloadSampleSet()
{
    if (sampleSet != NULL)
    {
        sampleSet->load(std::stoi((*samplesToSave).c_str()));
    }
}

void App::generateSampleSetPositions()
{
    // generate n random points and store their colors, compare to the color in the corresponding texture , add to a cost function
    int n = atoi(samplesToSave->c_str());
    Vector2 UV, IJ;
    Vector3 N, P;
    int selectedModel;

    for (int i = 0; i < n; i++)
    {
        clearAllActors();
        SceneSample ss = generateSceneSample(&selectedModel, &P, &UV, &N, i);
        sampleSet->addSample(ss);

    }
    sampleSet->save();
}

void App::generateSampleSetValues()
{
    std::stringstream args;
    runPythonScriptFromDataFiles("IrradianceSensorRender.py", std::string(m_scene->m_name.c_str()) + " " + sampleSet->m_sampleSetName + " " + std::string((*samplesToSave).c_str()), true, true);
    sampleSet->load(std::stoi((*samplesToSave).c_str()));
}

void App::generateSampleSetValuesFromProbes()
{
    int numCoeffs = std::atoi(optimizationSHBand.c_str());
    int numSamples = std::stoi((*samplesToSave).c_str());
    sampleSet->generateRGBValuesFromProbes(numSamples, numCoeffs);
}

void App::switchEditProbeStructure()
{

    if (!bManipulateProbesEnabled)
    {
		if (m_probeStructure->eType() == EProbeStructureType::Trilinear)
		{
			addWidget(m_probeStructure->getProbe(0)->getManipulator());

		}
		else
		{
			for (Probe* p : m_probeStructure->probeList)
			{
				addWidget(p->getManipulator());
			}
		}

    }
    else
    {
		if (m_probeStructure->eType() == EProbeStructureType::Trilinear)
		{
			removeWidget(m_probeStructure->getProbe(0)->getManipulator());
		}
		else
		{
			for (Probe* p : m_probeStructure->probeList)
			{
				removeWidget(p->getManipulator());
			}
		}
    }

	bManipulateProbesEnabled = !bManipulateProbesEnabled;
}

void App::updateProbeStructure()
{
	if (probeStructureLoaded() && bManipulateProbesEnabled)
	{
		switchEditProbeStructure();
	}

    G3D::String selectedScene = selectedSceneName();
    G3D::String probeStructureName = scenePane.probeStructureList->selectedValue();
    initializeProbeStructure(selectedScene, probeStructureName);
    bShouldUpdateProbeStructurePane = true;
}

void App::loadPreviousProbeStructure()
{
    if (previousProbeStructure != "")
    {
        initializeProbeStructure(selectedSceneName(), previousProbeStructure);
    }
    if (sampleSet)
    {
        sampleSet->probeStructure = m_probeStructure;
    }
}

void App::updateSampleSet()
{
    G3D::String selectedScene = selectedSceneName();
    G3D::String sampleSetName = scenePane.sampleSetList->selectedValue();

    //loadSurfaceSamples(selectedScene, sampleSetName);
    int numSamples = std::stoi((*samplesToSave).c_str());
    if (sampleSet != NULL)
    {
        delete sampleSet;
    }
	sampleSet = new SceneSampleSet(std::string(selectedScene.c_str()), std::string(sampleSetName.c_str()), m_scene->m_scale, numSamples);
	if (m_probeStructure != NULL)
	{
		sampleSet->probeStructure = m_probeStructure;
	}
}

void App::onAI()
{
    GApp::onAI();

    if (bShouldUpdateProbeStructurePane)
    {
        updateProbeStructurePane();
        bShouldUpdateProbeStructurePane = false;
    }

	if (actors.size() > 0)
	{
		screenPrintf("Actor pos: %s", actors[0].getPosition().toString().c_str());
	}

	if (numPassesLeft > 0)
	{
		tryOptimization();
        numPassesLeft--;

        if (numPassesLeft == 0)
        {
            popNotification("Optimization complete", "Finished all job!", 15);
        }
	}

    if (optimizing)
    {
        if (shouldAddAProbe)
        {
            findBestInitialConditions();
            shouldAddAProbe = false;
        }
        else
        {
            shouldAddAProbe = !tryOptimization();
        }

        if (shouldAddAProbe && (m_probeStructure->probeCount() == 6))
        {
            optimizing = false;
        }
    }

	if (m_probeStructure)
	{
		screenPrintf("ProbeStructure: %s", m_probeStructure->m_name.c_str());
	}
	
	//TODO: probe ray casting one day? if i find a real use for it
    for (int i = 0 ; i < actors.size(); i++)
	{
		if (interpolateCoefficients)
		{
			bool keepProbes = (i == 0);
			//actors[i].coefficients = triLinearInterpolation(i, &probeIndices, &weights, keepProbes);

			ProbeInterpolationRecord iRec = m_probeStructure->getInterpolationProbeIndicesAndWeights(actors[i].getPosition());

			//if (tetrahedralInterpolation(actors[i], &probeIndices, &weights))
			{
				actors[i].coefficients = getInterpolatedCoeffs(iRec, actors[i].drawBand);
			}
			
			if (showInterpolationProbes)
			{
				probesToRender.clear();
				for (int i = 0; i < iRec.probeIndices.size(); ++i)
				{
					probesToRender.push_back(m_probeStructure->getProbe(iRec.probeIndices[i]));
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


void App::offlineRender()
{

	std::stringstream command;
	command << "zcbox 1Probe Probes ";

	command << m_activeCamera->frame().translation.x << " ";
	command << m_activeCamera->frame().translation.y << " ";
	command << m_activeCamera->frame().translation.z << " ";
	command << m_activeCamera->frame().lookVector().x << " ";
	command << m_activeCamera->frame().lookVector().y << " ";
	command << m_activeCamera->frame().lookVector().z << " ";
	command << actors[0].getPosition().x << " ";
	command << actors[0].getPosition().y << " ";
	command << actors[0].getPosition().z << " ";
	command << offlineRenderingOptions.numSamples.c_str() << " ";
	command << offlineRenderingOptions.width.c_str() << " ";
	command << offlineRenderingOptions.height.c_str() << " ";
	command << integratorList[offlineRenderingOptions.integratorIndex].c_str() << " ";
	command << offlineRenderingOptions.gamma.c_str() << " ";
	command << filmTypeList[offlineRenderingOptions.filmTypeIndex].c_str() << " ";

	debugPrintf("%s\n", command.str().c_str());

	runPythonScriptFromDataFiles("CreateSceneRender.py", command.str(), true, true);
}

/*
//	GUI
*/

G3D::Array<G3D::String> getFoldersInFolder(const G3D::String& path)
{
	G3D::Array<G3D::String> folderList;
	FileSystem::ListSettings ls;
	ls.includeParentPath = false;
	ls.recursive = false;
	FileSystem::clearCache("");

	FileSystem::list(path + "/*", folderList, ls);

	return folderList;
}

void App::copyProbeStructure(String& sceneName, String& sourceProbeStructureName, String& dstProbeStructureName)
{
    String srcPath = probeStructureFoldersPath() + "/" + sourceProbeStructureName;
    String dstPath = probeStructureFoldersPath() + "/" + dstProbeStructureName;

    copyDir(srcPath, dstPath);
}

void App::createNewProbeStructure(String& sceneName, String& newProbeStructureName)
{
    String probeStructurePath = "../data-files/Scenes/" + sceneName + "/ProbeStructures/" + newProbeStructureName;

    createFolder(probeStructurePath);
    createFolder(probeStructurePath + "/Probes");
    createFolder(probeStructurePath + "/Normals");
    createFolder(probeStructurePath + "/Positions");
    createFolder(probeStructurePath + "/Coefficients");
    createEmptyFile(probeStructurePath + "/info.txt");

    std::fstream infoFile( (probeStructurePath + "/info.txt").c_str() , std::fstream::out);
    infoFile << "type wNN"     << std::endl;
    infoFile << "dimensions 0" << std::endl;
    infoFile << "width 128"    << std::endl;
    infoFile << "height 64"    << std::endl;
    infoFile << "gamma 1"      << std::endl;
    infoFile << "type wNN"     << std::endl;

    createEmptyFile(probeStructurePath + "/probeList.txt");
}

void App::createCopyProbeStructureWindow()
{
    windowCopyProbeStructure = GuiWindow::create("Copy probe structure");
    GuiPane* pane = windowCopyProbeStructure->pane();

    pane->beginRow();
    pane->addTextBox("Name: ", &newProbeStructureOptions.name);
    pane->endRow();

    pane->addButton("Ok", [this]()
    {
        copyProbeStructure(m_scene->m_name, m_probeStructure->m_name, newProbeStructureOptions.name);
        newProbeStructureOptions.name = "";
        windowCopyProbeStructure->setVisible(false);
    });
    pane->addButton("Cancel", [this]()
    {
        newProbeStructureOptions.name = "";
        windowCopyProbeStructure->setVisible(false);
    });
    windowCopyProbeStructure->pack();
    addWidget(windowCopyProbeStructure);
    windowCopyProbeStructure->moveTo(Point2(500, 500));
    windowCopyProbeStructure->setVisible(false);
}

void App::createNewProbeStructureWindow()
{
	windowNewProbeStructure = GuiWindow::create("New probe structure");
	GuiPane* pane = windowNewProbeStructure->pane();
    
    pane->beginRow();
    pane->addTextBox("Name: ", &newProbeStructureOptions.name);
    pane->endRow();

	pane->addButton("Ok", [this]()
	{
		createNewProbeStructure(m_scene->m_name, newProbeStructureOptions.name);
        newProbeStructureOptions.name = "";
        windowNewProbeStructure->setVisible(false);
    });
	pane->addButton("Cancel", [this]()
	{
        newProbeStructureOptions.name = "";
		windowNewProbeStructure->setVisible(false);
	});
	windowNewProbeStructure->pack();
	addWidget(windowNewProbeStructure);
	windowNewProbeStructure->moveTo(Point2(500, 500));
	windowNewProbeStructure->setVisible(false);
}

void App::createNewSampleSet(String& sceneName, String& newSampleSetName)
{
	G3D::String newSSPath = "../data-files/Scenes/" + sceneName + "/SampleSets/" + newSampleSetName;
	createFolder(newSSPath.c_str());
	createFolder((newSSPath + "/irradiance").c_str());
	createEmptyFile((newSSPath + "/SamplePositions.txt").c_str());
	createEmptyFile((newSSPath + "/IrradianceResults2.txt").c_str());
}


void App::createNewProbeWindow()
{
	windowNewProbe = GuiWindow::create("New probe");

	GuiPane* pane = windowNewProbe->pane();
	pane->addTextBox("Position: ", &sNewProbePosition);

	pane->addButton("Ok", [this]()
	{
		Vector3 newPosition = StringToVector3(sNewProbePosition);

		m_probeStructure->addProbe(newPosition);
		windowNewProbe->setVisible(false);
        bShouldUpdateProbeStructurePane = true;

	});

	pane->addButton("Cancel", [this]()
	{
		sNewProbePosition = "";
		windowNewProbe->setVisible(false);
	});

	windowNewProbe->pack();
	addWidget(windowNewProbe);
	windowNewProbe->moveTo(Point2(500, 500));
	windowNewProbe->setVisible(false);
}

void App::createRenameOptimizationWindow()
{
    windowRenameOptimization = GuiWindow::create("New sample set");
    GuiPane* pane = windowRenameOptimization->pane();
    pane->addTextBox("Name: ", &sCurrentOptimizationName);
    pane->addButton("Ok", [this]()
    {
        String currentOptimizationPath = currentOptimizationFolderPath();

        String newName = g3dString(m_currentOptimization) + " - " + sCurrentOptimizationName;
        String newOptimizationPath = optimizationFolderPath() + "/" + newName;
        std::rename(currentOptimizationPath.c_str(), newOptimizationPath.c_str());
    });
    pane->addButton("Cancel", [this]()
    {
        windowRenameOptimization->setVisible(false);
    });
    windowRenameOptimization->pack();
    addWidget(windowRenameOptimization);
    windowRenameOptimization->moveTo(Point2(500, 500));
    windowRenameOptimization->setVisible(false);
}


void App::createNewSampleSetWindow()
{
	windowNewSampleSet = GuiWindow::create("New sample set");
	GuiPane* pane = windowNewSampleSet->pane();
	pane->addTextBox("Name: ", &sNewSampleSetName);
	pane->addButton("Ok", [this]()
	{
		createNewSampleSet(m_scene->m_name, sNewSampleSetName);
		windowNewSampleSet->setVisible(false);
		generateSampleSetList();
		sNewSampleSetName = "";
	});
	pane->addButton("Cancel", [this]()
	{
		sNewSampleSetName = "";
		windowNewSampleSet->setVisible(false);
	});
	windowNewSampleSet->pack();
	addWidget(windowNewSampleSet);
	windowNewSampleSet->moveTo(Point2(500, 500));
	windowNewSampleSet->setVisible(false);
}

void App::createNewOptimizationSettings()
{
    String currentOptimizationPath = currentOptimizationFolderPath();
    const G3D::String sceneName = selectedSceneName();
    const G3D::String basePath = optimizationFolderPath() + "/";
    FileSystem::clearCache("");
    const String folderName(generateFolderNameBaseAnySuffix(basePath));

    m_currentOptimization++;

    createFolder(folderName.c_str());
    createFolder((folderName + "/screens").c_str());
    createEmptyFile((folderName + "/samples.txt").c_str());
    createEmptyFile((folderName + "/values.txt").c_str());

    //if (bKeepRefValuesOnNewOptimization)
    //{
    //    copyFile(currentOptimizationPath + "/ref_values.txt", currentOptimizationPath + "/ref_values.txt");
    //}
    //else
    {
        createEmptyFile((folderName + "/ref_values.txt").c_str());
    }
    createEmptyFile((folderName + "/errorlog.txt").c_str());
    createEmptyFile((folderName + "/dplog.txt").c_str());
    createEmptyFile((folderName + "/triplets.txt").c_str());
    createEmptyFile((folderName + "/log.txt").c_str());
}

void App::makeGui() {
	shared_ptr<GuiWindow> gui = GuiWindow::create("Parameters");
	GuiPane* pane = gui->pane();

	GuiTabPane* tabPane = pane->addTabPane();
	GuiPane* tab = tabPane->addTab("General Controls");

	tab->beginRow();
	tab->addCheckBox("Interpolate Coefficients", &interpolateCoefficients);
    tab->addButton("displace", GuiControl::Callback(this, &App::displaceProbes), GuiTheme::TOOL_BUTTON_STYLE);
    tab->addButton("New optimization", [this]()
    {
        createNewOptimizationSettings();
    }
    , GuiTheme::TOOL_BUTTON_STYLE);
    tab->addCheckBox("Keep ref values", &bKeepRefValuesOnNewOptimization);
    tab->addButton("Rename current optimization", [this]()
    {
        windowRenameOptimization->setVisible(true);
    }
    , GuiTheme::TOOL_BUTTON_STYLE);
	tab->addTextBox("Max probe step", &maxProbeStepLength);
	tab->addTextBox("SHBand", &optimizationSHBand);
	tab->addCheckBox("Update probes", &bUpdateProbesOnOptimizationPass);
	tab->endRow();

	tab->beginRow();
    tab->addButton("Compute samplesRGB", [this]() 
    {
        int numSamples = std::atoi(numOptimizationSamples.c_str());
        int numCoeffs = std::atoi(optimizationSHBand.c_str());
        String outputFile = currentOptimizationFolderPath() + "/values.txt";
        sampleSet->generateRGBValuesFromProbes(numSamples, numCoeffs, outputFile, 0);
    }
    , GuiTheme::TOOL_BUTTON_STYLE);

	tab->addButton("Compute ref_values", [this]() 
    {
        int numSamples = std::atoi(numOptimizationSamples.c_str());
        int numCoeffs = std::atoi(optimizationSHBand.c_str());
        String outputFile = currentOptimizationFolderPath() + "/ref_values.txt";
        if (bOptimizeWithMitsubaSamples)
        {
            sampleSet->generateRGBValuesFromSamples(numSamples, outputFile, 0);
        }
        else
        {
            sampleSet->generateRGBValuesFromProbes(numSamples, numCoeffs, outputFile, 0);
        }
        bTakeRefScreenshot = true;
    }
    ,GuiTheme::TOOL_BUTTON_STYLE);

    tab->addButton("Compute triplets", GuiControl::Callback(this, &App::computeTriplets), GuiTheme::TOOL_BUTTON_STYLE);

    tab->addButton("A", [this]() 
	{
		int numSamples = std::atoi(numOptimizationSamples.c_str());
        int numCoeffs = std::atoi(optimizationSHBand.c_str());
        sampleSet->outputWeightsMatrixToFile(numSamples, numCoeffs, currentOptimizationFolderPath());
    }
    , GuiTheme::TOOL_BUTTON_STYLE);

    tab->addButton("b", [this]() 
	{
		int numSamples = std::atoi(numOptimizationSamples.c_str());
        int numCoeffs = std::atoi(optimizationSHBand.c_str());
        sampleSet->outputBVectorToFile(numSamples, numCoeffs, currentOptimizationFolderPath());
    }
    , GuiTheme::TOOL_BUTTON_STYLE);

	tab->addButton("A*x", [this]()
	{
		int numSamples = std::atoi(numOptimizationSamples.c_str());
        int numCoeffs = std::atoi(optimizationSHBand.c_str());

		WeightMatrixType A = sampleSet->generateWeightsMatrix(numSamples, numCoeffs);

		String xPath = currentOptimizationFolderPath() + "/x.txt";
		std::fstream xFile(xPath.c_str(), std::fstream::in);
		std::string line;

		int numProbes = m_probeStructure->probeCount();
		Eigen::VectorXd x(numProbes * 3);
		int i = 0;
		while (std::getline(xFile, line))
		{
			x(i) = std::stof(line);
			i++;
		}
		xFile.close();

		Eigen::VectorXd b = A * x;

		String outputPath = currentOptimizationFolderPath() + "/AtimesX.txt";
		std::fstream output(outputPath.c_str(), std::fstream::out);
		output << b;
		output.close();

	}
	, GuiTheme::TOOL_BUTTON_STYLE);

	tab->addCheckBox("Path traced samples", &bOptimizeWithMitsubaSamples);
	tab->addButton("Optimization pass (GO!)", GuiControl::Callback(this, &App::startOptimizationPasses), GuiTheme::TOOL_BUTTON_STYLE);
    tab->addTextBox("Num passes", &tbNumPassesLeft);
    tab->addTextBox("Num samples", &numOptimizationSamples);
	tab->addCheckBox("log", &logSampleSet);
	tab->endRow();

	tab->beginRow();
	tab->addButton("Find Initial Conditions", GuiControl::Callback(this, &App::findBestInitialConditions), GuiTheme::TOOL_BUTTON_STYLE);
	tab->addTextBox("NumTries", &m_sNumICTries);
	tab->addTextBox("NumProbes", &m_sNumICProbes);
	tab->addButton("tryTotal", [this]()
	{
		
		findBestInitialConditions();
		shouldAddAProbe = false;
        optimizing = true;


	}, GuiTheme::TOOL_BUTTON_STYLE);
	tab->endRow();

	tab = tabPane->addTab("Display");
	tab->beginRow();
	tab->addCheckBox("Direct", &(bRenderDirect));
	tab->addCheckBox("Indirect (probes)", &(bRenderIndirect));
	tab->addCheckBox("* BRDF", &(bRenderMultiplyIndirectByBRDF));
	tab->addSlider("Mutiplier", &shadingMultiplier, 0.0f, 5.0f);
	tab->addCheckBox("Show Interp Probes", &showInterpolationProbes);
	tab->addCheckBox("Show ALL Probes", &showAllProbes);
	tab->addCheckBox("CPUInterpolation", &CPUInterpolation);
	tab->addCheckBox("highlight probes", &highlightProbes);
	tab->addTextBox("SH band (-1 = all)", &shadingSHBand);

	tab->endRow();

	tab->beginRow();


	tab->addCheckBox("AOF", &(bRenderAO));
	tab->addCheckBox("Shadow Maps", &(bRenderShadowMaps));
	tab->addTextBox("NumSamples", &maxSamplesPointsToDraw);
	tab->addCheckBox("Show", &showSamples);
	tab->addCheckBox("Show normals", &showSampleNormals);
	tab->addCheckBox("Show dark samples", &showDarkSamples);
	tab->addSlider("F", &sampleMultiplier, 1.0f, 50.f);
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
	for (int i = 0; i < actors.size(); i++)
	{
		tab->beginRow();
		tab->addLabel(GuiText(actors[i].name));
		tab->addCheckBox("visible", &actors[i].updateVisibility);
		tab->addSlider("Phong", &actors[i].phongExponent, 1.0f, 60.f);
		tab->addSlider("Max band", &actors[i].drawBand, 0, maxDrawBand);
		tab->addCheckBox("show this band only", &actors[i].showOnlyOneBand);
		tab->endRow();
	}

	tab = tabPane->addTab("Offline Render");

	tab->beginRow();
	tab->addDropDownList("Integrator", integratorList, &(offlineRenderingOptions.integratorIndex), [this]() {});
	tab->addDropDownList("Film Type", filmTypeList, &(offlineRenderingOptions.filmTypeIndex), [this]() {});
	tab->addTextBox("Gamma", &(offlineRenderingOptions.gamma));
	tab->addTextBox("NumSamples", &(offlineRenderingOptions.numSamples));
	tab->addTextBox("Width", &(offlineRenderingOptions.width));
	tab->addTextBox("Height", &(offlineRenderingOptions.height));
	tab->endRow();

	tab->beginRow();
	tab->addCheckBox("Show CMD", &(offlineRenderingOptions.showWindow));
	tab->addCheckBox("Requireclose CMD", &(offlineRenderingOptions.requireToCloseWindow));
	tab->addButton("OfflineRender", GuiControl::Callback(this, &App::offlineRender), GuiTheme::TOOL_BUTTON_STYLE);
	tab->endRow();



	addSampleSetPane(tabPane);

	probeStructurePane = tabPane->addTab("Probes");

	gui->pack();
	addWidget(gui);
	gui->moveTo(Point2(10, 50));
}

void App::generateSampleSetList()
{
	G3D::String selectedScene = selectedSceneName();
	G3D::Array<G3D::String> sampleSetList = getFoldersInFolder(sampleSetFoldersPath());

	scenePane.sampleSetList->setList(sampleSetList);
}

void App::updateProbeStructurePane()
{
	probeStructurePane->removeAllChildren();

	G3D::String selectedScene = selectedSceneName();

	probeStructurePane->beginRow();

	G3D::Array<G3D::String> probeStructureList = getFoldersInFolder("../data-files/Scenes/" + selectedScene + "/ProbeStructures");
	scenePane.probeStructureList = probeStructurePane->addDropDownList("ProbeStructure", probeStructureList, NULL, GuiControl::Callback(this, &App::updateProbeStructure));

    probeStructurePane->addButton(GuiText("New"), [this]() 
    {
        windowNewProbeStructure->setVisible(true);
        bShouldUpdateProbeStructurePane = true;
    }
    , GuiTheme::TOOL_BUTTON_STYLE);

    probeStructurePane->addButton(GuiText("New from copy..."), [this]()
    {
        if (!probeStructureLoaded())
        {
            return;
        }

        windowCopyProbeStructure->setVisible(true);
    }
    , GuiTheme::TOOL_BUTTON_STYLE);

    probeStructurePane->addButton(GuiText("Switch Back"), GuiControl::Callback(this, &App::loadPreviousProbeStructure), GuiTheme::TOOL_BUTTON_STYLE);
	probeStructurePane->addButton(GuiText("Edit mode"), GuiControl::Callback(this, &App::switchEditProbeStructure), GuiTheme::TOOL_BUTTON_STYLE);
	probeStructurePane->addButton(GuiText("Add probe"), [this]()
	{
		windowNewProbe->setVisible(true);
	}
	, GuiTheme::TOOL_BUTTON_STYLE);

	probeStructurePane->addButton(GuiText("Update all"), [this]()
	{
		m_probeStructure->savePositions(true);
		m_probeStructure->generateProbes("all");
		m_probeStructure->extractSHCoeffs();

		popNotification("Job finished", "Probe structure update complete", 15);
	}
	, GuiTheme::TOOL_BUTTON_STYLE);

    probeStructurePane->addButton(GuiText("Save positions"), [this]() 
    { 	
        m_probeStructure->savePositions(true);
    } 
    , GuiTheme::TOOL_BUTTON_STYLE);

    probeStructurePane->addButton(GuiText("Update Probes (all)"), [this]() { m_probeStructure->generateProbes("all"); }, GuiTheme::TOOL_BUTTON_STYLE);
    probeStructurePane->addButton(GuiText("Update Probes (probes)"), [this]() { m_probeStructure->generateProbes("Probes"); }, GuiTheme::TOOL_BUTTON_STYLE);
    probeStructurePane->addButton(GuiText("Update Probes (positions)"), [this]() { m_probeStructure->generateProbes("Positions"); }, GuiTheme::TOOL_BUTTON_STYLE);
    probeStructurePane->addButton(GuiText("Update Probes (normals)"), [this]() { m_probeStructure->generateProbes("Normals"); }, GuiTheme::TOOL_BUTTON_STYLE);
	probeStructurePane->addButton(GuiText("Extract coeffs"), [this]() { m_probeStructure->extractSHCoeffs(); }, GuiTheme::TOOL_BUTTON_STYLE);
    probeStructurePane->addCheckBox("Use gradients", &useSHGradients);
	probeStructurePane->endRow();

	probeStructurePane->beginRow();
	

	probeStructurePane->endRow();

	if (m_probeStructure != NULL)
	{
		int currentPSIndex = probeStructureList.findIndex(m_probeStructure->name());
		scenePane.probeStructureList->setSelectedIndex(currentPSIndex);


		probeStructurePane->beginRow();

		const G3D::String name = "Currently loaded : " + m_probeStructure->name();
		probeStructurePane->addLabel(name);
		const G3D::String type = "Type : " + m_probeStructure->type();
		probeStructurePane->addLabel(type);
		probeStructurePane->addLabel("Number of probes : " + String(std::to_string(m_probeStructure->probeCount())));
        probeStructurePane->endRow();

        probeStructurePanelOptions.gamma      = g3dString(m_probeStructure->gamma()     );
        probeStructurePanelOptions.numSamples = g3dString(m_probeStructure->numSamples());
        probeStructurePanelOptions.width      = g3dString(m_probeStructure->width()     );
        probeStructurePanelOptions.height     = g3dString(m_probeStructure->height()    );
        probeStructurePanelOptions.step     =   g3dString(m_probeStructure->m_step	    );
        
        probeStructurePane->beginRow();
        probeStructurePane->addTextBox("Gamma",      &(probeStructurePanelOptions.gamma));
        probeStructurePane->addTextBox("Width",      &(probeStructurePanelOptions.width));
        probeStructurePane->addTextBox("Height",     &(probeStructurePanelOptions.height));
		probeStructurePane->addTextBox("NumSamples", &(probeStructurePanelOptions.numSamples));
		probeStructurePane->addTextBox("Step",	     &(probeStructurePanelOptions.step));

		int currentIntegratorIndex = integratorList.findIndex(m_probeStructure->m_integrator);
		G3D::GuiDropDownList* list = probeStructurePane->addDropDownList("Integrator", integratorList, &(probeStructurePanelOptions.integratorIndex), [this]() 
		{ 
			String selected = integratorList[probeStructurePanelOptions.integratorIndex];
			m_probeStructure->setIntegrator(selected);
		});
		list->setSelectedIndex(currentIntegratorIndex);

		probeStructurePane->addButton(GuiText("Save infos"), [this]() 
		{
            m_probeStructure->setGamma(     std::stof(probeStructurePanelOptions.gamma.c_str()));
            m_probeStructure->setWidth(     std::stoi(probeStructurePanelOptions.width.c_str()));
            m_probeStructure->setHeight(    std::stoi(probeStructurePanelOptions.height.c_str()));
            m_probeStructure->setNumSamples(std::stoi(probeStructurePanelOptions.numSamples.c_str()));
			m_probeStructure->saveInfoFile(); 
		}, GuiTheme::TOOL_BUTTON_STYLE);
		probeStructurePane->addButton(GuiText("Delete all probes"), [this]()
		{
			m_probeStructure->deleteAllProbes();
		}, GuiTheme::TOOL_BUTTON_STYLE);

		probeStructurePane->endRow();
	}
	else
	{
		probeStructurePane->addLabel("No probe structure loaded");
	}

}


void App::addSampleSetPane(GuiTabPane* tabPane)
{
	GuiPane* tab = tabPane->addTab("Scene");

	G3D::Array<G3D::String> sceneList = getFoldersInFolder(scenesPath());

	scenePane.selectedSceneList = tab->addDropDownList("Scene", sceneList, NULL, GuiControl::Callback(this, &App::updateSelectedScenePane));

	if (sceneLoaded())
	{
		scenePane.selectedSceneList->setSelectedValue(m_scene->m_name);
	}

	G3D::String selectedScene = selectedSceneName();

	tab->beginRow();
	G3D::Array<G3D::String> sampleSetList = getFoldersInFolder(sampleSetFoldersPath());

	scenePane.sampleSetList = tab->addDropDownList("SampleSet", sampleSetList, NULL, GuiControl::Callback(this, &App::updateSampleSet));

	if (sampleSetLoaded())
	{
		scenePane.sampleSetList->setSelectedValue(String(sampleSet->m_sampleSetName.c_str()));
	}

	tab->addButton(GuiText("New"), [this]() { windowNewSampleSet->setVisible(true); }, GuiTheme::TOOL_BUTTON_STYLE);
	tab->addButton(GuiText("Clear values"), GuiControl::Callback(this, &App::clearSampleSetValues), GuiTheme::TOOL_BUTTON_STYLE);
	tab->addButton(GuiText("Clear positions"), GuiControl::Callback(this, &App::clearSampleSetPositions), GuiTheme::TOOL_BUTTON_STYLE);
	tab->addButton(GuiText("Reload"), GuiControl::Callback(this, &App::reloadSampleSet), GuiTheme::TOOL_BUTTON_STYLE);
	tab->addTextBox("samplesToSave", samplesToSave);
	tab->addButton("Generate Positions", GuiControl::Callback(this, &App::generateSampleSetPositions), GuiTheme::TOOL_BUTTON_STYLE);
    tab->addButton("Generate Values (mitsuba)", GuiControl::Callback(this, &App::generateSampleSetValues), GuiTheme::TOOL_BUTTON_STYLE);
    tab->addButton("Generate Values (probes)", GuiControl::Callback(this, &App::generateSampleSetValuesFromProbes), GuiTheme::TOOL_BUTTON_STYLE);
    tab->addCheckBox("write samples", &saveSample);
	tab->endRow();
}

void App::updateSelectedScenePane()
{
    G3D::String selectedScene = selectedSceneName();
	debugPrintf("update ldofgio");

	loadScene(selectedScene);

	G3D::Array<G3D::String> probeStructureList = getFoldersInFolder(probeStructureFoldersPath());
	scenePane.probeStructureList->setList(probeStructureList);

	G3D::Array<G3D::String> sampleSetList = getFoldersInFolder(sampleSetFoldersPath());
	scenePane.sampleSetList->setList(sampleSetList);
}

void App::onUserInput(UserInput* userInput)
{
    GApp::onUserInput(userInput);

    if (userInput->keyDown(GKey('u')))
    {
        bShouldUpdateProbeStructurePane = true;
    }
    else if (userInput->keyDown(GKey('0')))
    {
        numPassesLeft = 0;
    }
    else if (userInput->keyDown(GKey('p')))
    {
        showAllProbes = !showAllProbes;
    }
	
}

G3D::String App::scenesPath()
{
    return "../data-files/Scenes";
}

G3D::String App::selectedSceneName()
{
    return scenePane.selectedSceneList->selectedValue();
}

G3D::String App::optimizationFolderPath()
{
    return "../data-files/Scenes/" + selectedSceneName() + "/Optimizations";
}

G3D::String App::currentOptimizationFolderPath()
{
    return optimizationFolderPath() + "/" + g3dString(m_currentOptimization);
}

G3D::String App::probeStructureFoldersPath()
{
    return "../data-files/Scenes/" + selectedSceneName() + "/ProbeStructures";
}

G3D::String App::sampleSetFoldersPath()
{
    return "../data-files/Scenes/" + selectedSceneName() + "/SampleSets";
}

G3D::String App::loadedProbeStructurePath()
{
    String name = m_probeStructure->m_name;
    return probeStructureFoldersPath() + "/" + name;
}

bool App::probeStructureLoaded()
{
    return m_probeStructure != NULL;
}


bool App::sceneLoaded()
{
	return m_scene != NULL;
}

bool App::sampleSetLoaded()
{
	return sampleSet != NULL;
}

const char* optionFilePath()
{
	return "../data-files/options.txt";
}
void App::saveOptions()
{
	using json = nlohmann::json;
	json optionJSON;

	optionJSON["sceneName"] =						selectedSceneName().c_str();
	optionJSON["optimizationSHBand"] =				optimizationSHBand.c_str();
	optionJSON["numPassesLeft"] =					tbNumPassesLeft.c_str();
	optionJSON["numOptimizationSamples"] =			numOptimizationSamples.c_str();
	optionJSON["probeStructureName"] =				m_probeStructure->m_name.c_str();
	optionJSON["sampleSetName"] =					sampleSet->m_sampleSetName.c_str();

	optionJSON["bShowAllProbes"] =					showAllProbes;
	optionJSON["bUpdateProbesOnOptimizationPass"] = bUpdateProbesOnOptimizationPass;
	optionJSON["bRenderDirect"] =					bRenderDirect;
	optionJSON["bRenderIndirect"] =					bRenderIndirect;
	optionJSON["sampleMultiplier"] =				sampleMultiplier;

	std::fstream optionFile(optionFilePath(), std::fstream::out);
	optionFile << optionJSON;
	optionFile.close();
}

using json = nlohmann::json;

String App::loadStringOption(String name, nlohmann::json& optionJSON)
{
	try
	{
		json::string_t val = optionJSON[name.c_str()];
		return String(val.c_str());
	}
	catch (std::exception e)
	{
		return "";
	}
}

bool App::loadBoolOption(String name, nlohmann::json& optionJSON)
{
	try
	{
		json::boolean_t val = optionJSON[name.c_str()];
		return val;
	}
	catch (std::exception e)
	{
		return false;
	}
}

void App::loadOptions()
{
	json optionJSON;

	std::fstream optionFile(optionFilePath(), std::fstream::in);
	optionJSON << optionFile;
	optionFile.close();


	loadScene(loadStringOption("sceneName", optionJSON));
	bUpdateProbesOnOptimizationPass = loadBoolOption("bUpdateProbesOnOptimizationPass", optionJSON);
	showAllProbes = loadBoolOption("bShowAllProbes", optionJSON);
	bRenderDirect = loadBoolOption("bRenderDirect", optionJSON);
	bRenderIndirect = loadBoolOption("bRenderIndirect", optionJSON);

	try
	{
		json::string_t probeStructureName = optionJSON["probeStructureName"];
		initializeProbeStructure(m_scene->m_name, String(probeStructureName.c_str()));
		bShouldUpdateProbeStructurePane = true;
	}
	catch (std::exception e)
	{

	}

	try
	{
		json::string_t sampleSetName = optionJSON["sampleSetName"];
		sampleSet = new SceneSampleSet(std::string(m_scene->m_name.c_str()), std::string(sampleSetName.c_str()), m_scene->m_scale, -1);
		if (m_probeStructure != NULL)
		{
			sampleSet->probeStructure = m_probeStructure;
		}
	}
	catch (std::exception e)
	{

	}

	try
	{
		json::string_t val = optionJSON["optimizationSHBand"];
		optimizationSHBand = String(val.c_str());
	}
	catch (std::exception e)
	{

	}

	try
	{
		json::string_t val = optionJSON["numPassesLeft"];
		tbNumPassesLeft = String(val.c_str());
	}
	catch (std::exception e)
	{

	}

	try
	{
		json::string_t val = optionJSON["numOptimizationSamples"];
		numOptimizationSamples = String(val.c_str());
	}
	catch (std::exception e)
	{

	}
}

bool App::onEvent(const GEvent& event)
{
	if (event.type == GEventType::QUIT)
	{
		saveOptions();
		exit(0);
	}

	return false;
}