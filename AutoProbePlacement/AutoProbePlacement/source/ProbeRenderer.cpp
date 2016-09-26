#include "ProbeRenderer.h"

void ProbeRenderer::render
(RenderDevice*                       rd,
	const shared_ptr<Framebuffer>&      framebuffer,
	const shared_ptr<Framebuffer>&      depthPeelFramebuffer,
	LightingEnvironment&                lightingEnvironment,
	const shared_ptr<GBuffer>&          gbuffer,
	const Array<shared_ptr<Surface>>&   allSurfaces) {

	alwaysAssertM(!lightingEnvironment.ambientOcclusionSettings.enabled || notNull(lightingEnvironment.ambientOcclusion),
		"Ambient occlusion is enabled but no ambient occlusion object is bound to the lighting environment");

	const shared_ptr<Camera>& camera = gbuffer->camera();

	// Share the depth buffer with the forward-rendering pipeline
	framebuffer->set(Framebuffer::DEPTH, gbuffer->texture(GBuffer::Field::DEPTH_AND_STENCIL));
	if (notNull(depthPeelFramebuffer)) {
		depthPeelFramebuffer->resize(framebuffer->width(), framebuffer->height());
	}

	// Cull and sort
	Array<shared_ptr<Surface> > sortedVisibleSurfaces, forwardOpaqueSurfaces, forwardBlendedSurfaces;
	cullAndSort(gbuffer, allSurfaces, sortedVisibleSurfaces, forwardOpaqueSurfaces, forwardBlendedSurfaces);

	debugAssert(framebuffer);
	// Bind the main framebuffer
	rd->pushState(framebuffer); {
		rd->clear();
		rd->setProjectionAndCameraMatrix(camera->projection(), camera->frame());

		const bool needDepthPeel = lightingEnvironment.ambientOcclusionSettings.useDepthPeelBuffer && lightingEnvironment.ambientOcclusionSettings.enabled;
		computeGBuffer(rd, sortedVisibleSurfaces, gbuffer, needDepthPeel ? depthPeelFramebuffer : nullptr, lightingEnvironment.ambientOcclusionSettings.depthPeelSeparationHint);

		// Shadowing + AO
		computeShadowing(rd, allSurfaces, gbuffer, depthPeelFramebuffer, lightingEnvironment);

		// Maybe launch deferred pass
		if (deferredShading()) {
			renderProbeShading(rd, gbuffer, lightingEnvironment);
		}

		// Main forward pass
		renderOpaqueSamples(rd, deferredShading() ? forwardOpaqueSurfaces : sortedVisibleSurfaces, gbuffer, lightingEnvironment);

		// Prepare screen-space lighting for the *next* frame
		lightingEnvironment.copyScreenSpaceBuffers(framebuffer, gbuffer->colorGuardBandThickness(), gbuffer->depthGuardBandThickness());

		renderOpaqueScreenSpaceRefractingSamples(rd, deferredShading() ? forwardOpaqueSurfaces : sortedVisibleSurfaces, gbuffer, lightingEnvironment);

		// Samples that require blending
		if (m_orderIndependentTransparency) {
			renderOrderIndependentBlendedSamples(rd, forwardBlendedSurfaces, gbuffer, lightingEnvironment);
		}
		else {
			renderSortedBlendedSamples(rd, forwardBlendedSurfaces, gbuffer, lightingEnvironment);
		}

	} rd->popState();
}

void ProbeRenderer::renderProbeShading(RenderDevice* rd, const shared_ptr<GBuffer>& gbuffer, const LightingEnvironment& environment) {
	// Make a pass over the screen, performing shading
	rd->push2D(); {
		rd->setGuardBandClip2D(gbuffer->colorGuardBandThickness());

		// Don't shade the skybox on this pass because it will be forward rendered
		rd->setDepthTest(RenderDevice::DEPTH_GREATER);
		Args args;

		environment.setShaderArgs(args);
		gbuffer->setShaderArgsRead(args, "gbuffer_");

		args.setRect(rd->viewport());
		args.setMacro("myLolMacro", true);

		LAUNCH_SHADER("ProbeRenderer_deferredShade.pix", args);
	} rd->pop2D();
}