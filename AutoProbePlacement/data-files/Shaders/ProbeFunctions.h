#include "SH.h"



vec3 findProbe000(vec3 wsPos, vec3 normal, float step)
{
	vec3 toReturn;
	
	toReturn.x = floor( (wsPos.x + normal.x * EPSILON) / step ) * step;
	toReturn.y = floor( (wsPos.y + normal.y * EPSILON) / step ) * step;
	toReturn.z = floor( (wsPos.z + normal.z * EPSILON) / step ) * step;
	
	return toReturn;
}

uint findProbeIndex(vec3 probePosition, vec3 startingPosition, int dimensions[3], float step)
{
	vec3 fprobeSpaceCoords = (probePosition - startingPosition) / step;
	
	ivec3 probeSpaceCoords = ivec3( fprobeSpaceCoords );
	uint index = probeSpaceCoords.x * dimensions[1] * dimensions[2] +
				 probeSpaceCoords.y * dimensions[2] +
				 probeSpaceCoords.z;
	
	return index;
}

void findInterpolatingProbes(vec3 wsPos, vec3 normal, float step, out uint indices[8], out vec3 weights)
{
	vec3 probe000Pos = findProbe000(wsPos, normal, step);
	
	weights = (wsPos - probe000Pos) / step;
	
	indices[0] = findProbeIndex(probe000Pos, probeStructure.firstProbePosition, probeStructure.dimensions, step);
	
	ivec3 offsets = ivec3(probeStructure.dimensions[2] * probeStructure.dimensions[1], probeStructure.dimensions[2], 1);
	
	// TODO: fix the fact that this doesnt consider the borders.
	indices[1] = indices[0] + offsets[2]; 							 // p001
	indices[2] = indices[0] + offsets[1]; 							 // p010
	indices[3] = indices[0] + offsets[1] + offsets[2];               // p011
	indices[4] = indices[0] + offsets[0];                            // p100
	indices[5] = indices[0] + offsets[0] + offsets[2];               // p101
	indices[6] = indices[0] + offsets[0] + offsets[1];               // p110
	indices[7] = indices[0] + offsets[0] + offsets[1] + offsets[2];  // p111
}                                                           

vec3 getSHCoeffsOneProbe(uint probeIndex, uint shBand)
{
	return vec3( probeStructure.probes[probeIndex].coefficients[shBand * 3 + 0],
				 probeStructure.probes[probeIndex].coefficients[shBand * 3 + 1],
				 probeStructure.probes[probeIndex].coefficients[shBand * 3 + 2] );
}

vec3 trilinearInterpolation(vec3 values[8], vec3 weights)
{
	vec3 c00 = values[p000] * (1 - weights.x) + values[p100] * weights.x;
	vec3 c01 = values[p001] * (1 - weights.x) + values[p101] * weights.x;
	vec3 c10 = values[p010] * (1 - weights.x) + values[p110] * weights.x;
	vec3 c11 = values[p011] * (1 - weights.x) + values[p111] * weights.x;
	
	vec3 c0 = c00 * (1 - weights.y) + c10 * weights.y;
	vec3 c1 = c01 * (1 - weights.y) + c11 * weights.y;
	
	vec3 finalValue = c0 * (1 - weights.z) + c1 * weights.z;
	return finalValue;
}

vec3 interpolateCoefficientsOneBand(uint probeIndices[8], vec3 weights, uint shBand)
{
	vec3 values[8];
	
	for (int i = 0; i < 8; ++i)
	{
		values[i] = getSHCoeffsOneProbe(probeIndices[i], shBand);
	}
	
	return trilinearInterpolation(values, weights);
}

float getSHCoeff(float u, float v, uint shBand)
{
	vec4 value;
	switch (shBand)
	{
		case 0: value = texture( uSampler0, vec2(u,v) ) * 2 - vec4(1.0);
				break;
		case 1: value = texture( uSampler1, vec2(u,v) ) * 2 - vec4(1.0);
				break;
		case 2: value = texture( uSampler2, vec2(u,v) ) * 2 - vec4(1.0);
				break;
		case 3: value = texture( uSampler3, vec2(u,v) ) * 2 - vec4(1.0);
				break;
		case 4: value = texture( uSampler4, vec2(u,v) ) * 2 - vec4(1.0);
				break;
		case 5: value = texture( uSampler5, vec2(u,v) ) * 2 - vec4(1.0);
				break;
		case 6: value = texture( uSampler6, vec2(u,v) ) * 2 - vec4(1.0);
				break;
		case 7: value = texture( uSampler7, vec2(u,v) ) * 2 - vec4(1.0);
				break;
		case 8: value = texture( uSampler8, vec2(u,v) ) * 2 - vec4(1.0);
				break;
	}
	return value.x;
}


vec3 additiveValueOneBand(uint probeIndices[8], vec3 weights, uint shBand, float u, float v)
{
	vec3 interpolatedCoefficients = interpolateCoefficientsOneBand(probeIndices, weights, shBand);
	
	float shValue = getSHCoeff(u, v, shBand);
	if (shBand == 0) shValue = 0.2820947;
	
	return interpolatedCoefficients * shValue * phongCoeffs( kToLM(int(shBand))[0], r );
}

vec3 fullSHInterpolation(uint probeIndices[8], vec3 weights, float u, float v)
{
	vec3 toReturn = vec3(0, 0, 0);
	
	for (int shBand = 0; shBand < 9; ++shBand)
	{
		toReturn += additiveValueOneBand(probeIndices, weights, shBand, u, v);
	}
	return toReturn;
	
}

bool nearBorder(in vec3 interpolationWeights, in vec3 normal)
{
	bool nearBorder = false;

	vec3 threshold = vec3(0.01);

	if ( any( lessThan(interpolationWeights, threshold) ) && any( greaterThan(interpolationWeights, vec3(1) - threshold) ))
	{
		// if ( any( lessThan(normal, threshold) )  )
		{
			nearBorder = true;
		}
	}
	// if ( any( greaterThan(interpolationWeights, 1.f - threshold) ) )
	// {
		// nearBorder = true;
	// }
	return nearBorder;
}