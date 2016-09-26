#pragma once

#include "G3D\G3DAll.h"

class ProbeRenderer : public DefaultRenderer
{
public:
	static shared_ptr<Renderer> create() {
		return shared_ptr<ProbeRenderer>(new ProbeRenderer());
	}

	virtual void render
	(RenderDevice*                       rd,
		const shared_ptr<Framebuffer>&      framebuffer,
		const shared_ptr<Framebuffer>&      depthPeelFramebuffer,
		LightingEnvironment&                lightingEnvironment,
		const shared_ptr<GBuffer>&          gbuffer,
		const Array<shared_ptr<Surface>>&   allSurfaces) override;

	void renderProbeShading(RenderDevice* rd, const shared_ptr<GBuffer>& gbuffer, const LightingEnvironment& environment);
};

