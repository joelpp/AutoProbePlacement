
vec3 findNode000(vec3 wsPos, vec3 normal, float step)
{
	vec3 toReturn;
	
	toReturn.x = floor( wsPos.x / step ) * step;
	toReturn.y = floor( wsPos.y / step ) * step;
	toReturn.z = floor( wsPos.z / step ) * step;
	
	return toReturn;
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
