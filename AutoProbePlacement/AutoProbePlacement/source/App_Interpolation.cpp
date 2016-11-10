#include "App.h"

bool App::isUVInsideTriangle(Tri &t, Vector2 testUV, Vector3 &barycentric){
	const CPUVertexArray::Vertex &v0 = t.vertex(cva, 0);
	const CPUVertexArray::Vertex &v1 = t.vertex(cva, 1);
	const CPUVertexArray::Vertex &v2 = t.vertex(cva, 2);

	//debugPrintf("v0: %s\n",v0.texCoord0.toString().c_str());
	//debugPrintf("v1: %s\n",v1.texCoord0.toString().c_str());
	//debugPrintf("v2: %s\n",v2.texCoord0.toString().c_str());

	Vector3 u = Vector3(v0.texCoord0.x, v1.texCoord0.x, v2.texCoord0.x);
	Vector3 v = Vector3(v0.texCoord0.y, v1.texCoord0.y, v2.texCoord0.y);

	if (!pointInTriangleBB(testUV, u, v)) return false;

	barycentric = getUVBarycentricCoordinates(testUV, u, v);

	//debugPrintf("1Address of barycentric: %d\n", &barycentric);
	//debugPrintf("Weights: %s\n", barycentric.toString().c_str());

	return (!badWeights(barycentric,1.0));
}

bool App::pointInTriangleBB(Point2 UV, Vector3 u, Vector3 v){

		bool allOver = (UV.x > u.x) && (UV.x > u.y) && (UV.x > u.z);
		bool allUnder = (UV.x < u.x) && (UV.x < u.y) && (UV.x < u.z);
		bool allLeft = (UV.y < v.x) && (UV.y < v.y) && (UV.y < v.z);
		bool allRight = (UV.y > v.x) && (UV.y > v.y) && (UV.y > v.z);
		
		return !(allOver || allUnder || allLeft || allRight);	
}

Vector3 App::getUVBarycentricCoordinates(Point2 testUV, Vector3 u, Vector3 v){

    Matrix2 front = Matrix2(u.x - u.z, u.y - u.z, v.x - v.z, v.y-v.z);
    front = front.inverse();

    Vector2 back = Vector2(testUV.x - u.z, testUV.y - v.z);

    Vector2 ab = front * back;

	Vector3 toReturn = Vector3(ab.x, ab.y, 1 - ab.x - ab.y);
	//debugPrintf("0Address of toReturn: %d\n", &toReturn);
    return toReturn;
}
Vector3 App::triangleInterpolateVector3(Vector3 weights, Vector3 v0, Vector3 v1, Vector3 v2){

    return (weights.x * v0 + weights.y * v1 + weights.z * v2);
}

Array<Vector3> App::getInterpolatedCoeffs(ProbeInterpolationRecord iRec, int maxBand){
	Array<Vector3> interpolatedCoeffs = Array<Vector3>();
    for (int i = 0; i <= maxBand; i++){
        interpolatedCoeffs.append(Vector3(0,0,0));
    }

	for (int i = 0; i < iRec.probeIndices.size(); i++)
	{
		Probe* p = m_probeStructure->getProbe(iRec.probeIndices[i]);
		for (int j = 0; j <= maxBand; ++j)
		{
			interpolatedCoeffs[j] += p->coeffs[j] * iRec.weights[i];
		}
    }

	return interpolatedCoeffs;
}


Vector3 App::getBarycentricCoordinates(Point3 testPoint, Tetrahedron* t){
    Matrix3 part1;
	//debugPrintf("getting bar coords for tet %s", t->toString().c_str());
	Probe* probe0;
	Probe* probe1;
	Probe* probe2;
	Probe* probe3;

	if (actors[0].tetrahedronIndex == -1) return Vector3(-1,-1,-1);
	try{
		probe0 = m_probeStructure->getProbe(t->indexes[0]);
		probe1 = m_probeStructure->getProbe(t->indexes[1]);
		probe2 = m_probeStructure->getProbe(t->indexes[2]);
		probe3 = m_probeStructure->getProbe(t->indexes[3]);
	}
	catch(...){
		throw std::invalid_argument("you suck");
	}

    Vector3 v0 = probe0->position - probe3->position;
    Vector3 v1 = probe1->position - probe3->position;
    Vector3 v2 = probe2->position - probe3->position;
    Vector3 v3 = testPoint - probe3->position;

    part1 = Matrix3::fromColumns(v0, v1, v2);
    part1 = part1.inverse();
    Vector3 abc = part1 * v3;
    return abc;
}


bool App::tetrahedralInterpolation(Actor& actor, Array<int> *_probeIndices, Array<float> *_coeffs){

    Vector3 testPoint = actor.getPosition();

	Tetrahedron* t = m_probeStructure->getTetrahedron(actor.tetrahedronIndex);

        Vector3 weights;
        float d;
		float pDotFaceNormal = (testPoint - m_probeStructure->getFaceNormalPosition(actor.tetrahedronIndex)).dot(m_probeStructure->getFaceNormal(actor.tetrahedronIndex));
        // debugPrintf("p dot face normal = %f\n",pDotFaceNormal);
        
		// If the dot product is negative we are inside the hull
        if (pDotFaceNormal <= 0)
		{

			//remove particles for triangle visualisation
            //resetParticles();
            smallestBaryCoord = -1;

			// Calculate barycentric weights for current tetrahedron
            try
			{
				weights = getBarycentricCoordinates(testPoint, t);
			}
			catch(...)
			{
				throw std::invalid_argument("you suck!");
			}
            d = 1.0f - weights.sum();

			// If any of the weights are negative we are out of the current tetrahedron
            if (badWeights(weights, d))
			{
				// Switch the active tetrahedron and get the new weights
                t = updateNeighbor(weights, d, actor, t);

				// Get teh new weights
				try
				{
					weights = getBarycentricCoordinates(testPoint, t);
				}
				catch(...)
				{
					throw std::invalid_argument("you suck!");
				}
				d = 1.0f - weights.sum();
            }
			if (weights == Vector3(-1,-1,-1)) return false;
        }
        else
		{ // we are outside the hull
             screenPrintf("testing for tetrahedron = %s.\n", t->toString().c_str());

			 // find the T at which our testPoint is coplanar with the extruded triangle face for the current Tetrahedron
            double T = findTforCurrentTriangleFace(testPoint, t);
             screenPrintf("shitty t came back from asshole function as t = %f.\n", T);

			 // compute the locations of the extruded triangle's vertices
            Array<Vector3> triangleVertices = getTriangleVertices(T, t);

			// get the barycentric weights for our testpoint on the triangle face
            weights = findTriangleBarycentricCoordinates(testPoint, triangleVertices[0], triangleVertices[1], triangleVertices[2]);
			d = 1.0f - weights.x - weights.y - weights.z;
            screenPrintf("weights = %s.\n", weights.toString().c_str());
            // debugPrintf("Barycentric coords something like = %s.\n", barycentricCoordinates.toString().c_str());
            // debugPrintf("their sum = %f.\n", barycentricCoordinates.sum());
            
			// if any of the weights are negative we are extruding the wrong tetrahedron
            if (badWeights(weights, 1.0))
			{

				// so we have to update the active tetrahedron
                t = updateNeighbor(weights, 1.0, actor, t);

				// and find the T again for THAT triangle face
                T = findTforCurrentTriangleFace(testPoint, t);

				// and find the triangle vertex locations again
                triangleVertices = getTriangleVertices(T, t);

				// and find the weights again
                weights = findTriangleBarycentricCoordinates(testPoint, triangleVertices[0], triangleVertices[1], triangleVertices[2]);
                
            }
            smallestBaryCoord = minCoord(weights);

			// since we are outside the hull the 4th weight is 0
            d = 0;
        }

        _probeIndices->append(t->indexes[0]);
        _probeIndices->append(t->indexes[1]);
        _probeIndices->append(t->indexes[2]);
        _probeIndices->append(t->indexes[3]);
		
		_coeffs->append(weights.x);
		_coeffs->append(weights.y);
		_coeffs->append(weights.z);
		_coeffs->append(d);


		return true;
}

void App::tetrahedralInterpolation(G3D::Vector3 testPoint, Array<int> *probeIndices, Array<float> *weights)
{
	G3D::Vector3 spawnPoint = G3D::Vector3(2,2,2);

	Actor tempActor = Actor("temp", sphereModel, spawnPoint, shared_ptr<Texture>(), true);

	int numIterations = 10;
	G3D::Vector3 increment = (testPoint - spawnPoint) / 10.0f;

	for (int i = 0; i < numIterations; ++i)
	{
		tempActor.setPosition(tempActor.getPosition() + increment);
		tetrahedralInterpolation(tempActor, probeIndices, weights);
	}
}

int App::minCoord(Vector3 weights){
    if ((weights.x < weights.y) && (weights.x < weights.z)) return 0;
    else if (weights.y < weights.z) return 1;
    else return 2;
}

Array<Vector3> App::getTriangleVertices(double t, Tetrahedron* tet){

    //resetParticles();

	Array<int> indices = Array<int>();
	for (int i = 0; i < 4; i++)
		if (tet->neighbors[i] == -1) continue;
		else indices.append(tet->indexes[i]);

    Vector3 V0 = m_probeStructure->getProbe(indices[0])->position;
    Vector3 V1 = m_probeStructure->getProbe(indices[1])->position;
    Vector3 V2 = m_probeStructure->getProbe(indices[2])->position;
													
	Vector3 N0 = m_probeStructure->getProbe(indices[0])->normal;
    Vector3 N1 = m_probeStructure->getProbe(indices[1])->normal;
    Vector3 N2 = m_probeStructure->getProbe(indices[2])->normal;

    Array<Vector3> toReturn = Array<Vector3>();

    toReturn.append(V0 + t * N0);
    toReturn.append(V1 + t * N1);
    toReturn.append(V2 + t * N2);

    // debugPrintf("toReturn[0]: %s", toReturn[0].toString().c_str());

    // printVector3Array(toReturn);

	// For debugging/visualisation. Add spheres showing the vertices of the extrapolation triangle.
    //addParticle(toReturn[0].x, toReturn[0].y, toReturn[0].z);
    //addParticle(toReturn[1].x, toReturn[1].y, toReturn[1].z);
    //addParticle(toReturn[2].x, toReturn[2].y, toReturn[2].z);

	extrapolationTriangleVertices = toReturn;

    return toReturn;
}

Tetrahedron* App::updateNeighbor(Vector3 weights, float d, Actor& actor, Tetrahedron* t){
    int newIndex;
	float totalWeights[4] = { weights[0], weights[1], weights[2], d };
    //if (weights.x < 0){ newIndex = t->neighbors[0]; }
    //else if (weights.y < 0){ newIndex = t->neighbors[1]; }
    //else if (weights.z < 0){ newIndex = t->neighbors[2]; }
    //else if (d < 0){ newIndex = t->neighbors[3]; }
    //else{ }

	for (int i = 0 ; i < 4; ++i)
	{
		if (totalWeights[i] < 0)
		{
			newIndex = t->neighbors[i];
		
			if (newIndex == -1)
			{
				newIndex = t->neighbors[i+1];
			}
			break;
		}
	}


    actor.tetrahedronIndex = newIndex;

	if (newIndex == -1)
	{
		debugPrintf("Trying to access tetrahedron index -1.\n");
		debugPrintf("T: %s\n", t->toString().c_str());
		debugPrintf("Weights: %s\n", weights.toString().c_str());
		exit(1);
	}
    return m_probeStructure->getTetrahedron(newIndex);
}

bool App::badWeights(Vector3 v, float d){
	double threshold = -0.00;
    return ( ((v.x < threshold) || (v.y < threshold) || (v.z < threshold) || ( d < threshold )) ||
			((isNaN(v.x)) || (isNaN(v.y)) || (isNaN(v.z)) || (isNaN(d)) ));
}
    
double App::findTforCurrentTriangleFace(Point3 testPosition, Tetrahedron* tet){

	Array<int> indices = Array<int>();
	for (int i = 0; i < 4; i++)
		if (tet->neighbors[i] == -1) continue;
		else indices.append(tet->indexes[i]);

    Vector3 V0 = m_probeStructure->getProbe(indices[0])->position - m_probeStructure->getProbe(indices[2])->position;
    Vector3 V1 = m_probeStructure->getProbe(indices[1])->position - m_probeStructure->getProbe(indices[2])->position;
    Vector3 V2 = m_probeStructure->getProbe(indices[2])->position - testPosition;
												   
    Vector3 N0 = (m_probeStructure->getProbe(indices[0])->normal - m_probeStructure->getProbe(indices[2])->normal).fastUnit();
    Vector3 N1 = (m_probeStructure->getProbe(indices[1])->normal - m_probeStructure->getProbe(indices[2])->normal).fastUnit();
    Vector3 N2 = (m_probeStructure->getProbe(indices[2])->normal).fastUnit();

	screenPrintf("m_probeStructure->getProbe(%d).normal: %s\n", indices[0], m_probeStructure->getProbe(indices[0])->normal.toString().c_str());
	screenPrintf("m_probeStructure->getProbe(%d).normal: %s\n", indices[1], m_probeStructure->getProbe(indices[1])->normal.toString().c_str());
	screenPrintf("m_probeStructure->getProbe(%d).normal: %s\n", indices[2], m_probeStructure->getProbe(indices[2])->normal.toString().c_str());

    double t = extrapolationT;
    double prevT = extrapolationT + 100;

	int nIterations = 0;

    while ((fabs((t - prevT) / t) > 0.1) && nIterations < 100){
        prevT = t;
        t = t - extrapolationMatrixDeterminantFunction(t, V0, V1, V2, N0, N1, N2) / extrapolationMatrixDeterminantFunctionDerived(t, V0, V1, V2, N0, N1, N2);
		nIterations++;
	}
    extrapolationT = t;
    return t;
}

double App::extrapolationMatrixDeterminantFunction(double t, Vector3 V0, Vector3 V1, Vector3 V2, Vector3 N0, Vector3 N1, Vector3 N2){

    double a0 = N1.x * N2.y * N0.z + N0.x * N1.y * N2.z + N2.x * N0.y * N1.z - N2.x * N1.y * N0.z - N0.x * N2.y * N1.z - N1.x * N0.y * N2.z;

    double a1 =  V0.x * N1.y * N2.z + N0.x * V1.y * N2.z + N0.x * N1.y * V2.z + V1.x * N2.y * N0.z + N1.x * V2.y * N0.z + N1.x * N1.y * V0.z 
                + V2.x * N0.y * N1.z + N2.x * V0.y * N1.z + N2.x * N0.y * V1.z - V2.x * N1.y * N0.z - N2.x * V1.y * N0.z - N2.x * N1.y * V0.z
                 - V0.x * N2.y * N1.z - N0.x * V2.y * N1.z - N0.x * N2.y * V1.z - V1.x * N0.y * N2.z - N1.x * V0.y * N2.z - N1.x * N0.y * V2.z;   

    double a2 = V0.x * V1.y * N2.z + V0.x * N1.y * V2.z + N0.x * V1.y * V2.z + V1.x * V2.y * N0.z + V1.x * N2.y * V0.z + N1.x * V2.y * V0.z +
                V2.x * V0.y * N1.z + V2.x * N0.y * V1.z + N2.x * V0.y * V1.z - V2.x * V1.y * N0.z - V2.x * N1.y * V0.z - N2.x * V1.y * V0.z
                 - V0.x * V2.y * N1.z - V0.x * N2.y * V1.z - N0.x * V2.y * V1.z - V1.x * V0.y * N2.z - V1.x * N0.y * V2.z - N1.x * V0.y * V2.z;

    double a3 = V0.x * V1.y * V2.z + V1.x * V2.y * V0.z + V2.x * V0.y * V1.z - V2.x * V1.y * V0.z - V0.x * V2.y * V1.z - V1.x * V0.y * V2.z;

    return t*t*t * a0 + t*t * a1 + t * a2 + a3;
}


double App::extrapolationMatrixDeterminantFunctionDerived(double t, Vector3 V0, Vector3 V1, Vector3 V2, Vector3 N0, Vector3 N1, Vector3 N2)
{

    double a0 = N1.x * N2.y * N0.z + N0.x * N1.y * N2.z + N2.x * N0.y * N1.z - N2.x * N1.y * N0.z - N0.x * N2.y * N1.z - N1.x * N0.y * N2.z;

    double a1 =  V0.x * N1.y * N2.z + N0.x * V1.y * N2.z + N0.x * N1.y * V2.z + V1.x * N2.y * N0.z + N1.x * V2.y * N0.z + N1.x * N1.y * V0.z 
                + V2.x * N0.y * N1.z + N2.x * V0.y * N1.z + N2.x * N0.y * V1.z - V2.x * N1.y * N0.z - N2.x * V1.y * N0.z - N2.x * N1.y * V0.z
                 - V0.x * N2.y * N1.z - N0.x * V2.y * N1.z - N0.x * N2.y * V1.z - V1.x * N0.y * N2.z - N1.x * V0.y * N2.z - N1.x * N0.y * V2.z;   

    double a2 = V0.x * V1.y * N2.z + V0.x * N1.y * V2.z + N0.x * V1.y * V2.z + V1.x * V2.y * N0.z + V1.x * N2.y * V0.z + N1.x * V2.y * V0.z +
                V2.x * V0.y * N1.z + V2.x * N0.y * V1.z + N2.x * V0.y * V1.z - V2.x * V1.y * N0.z - V2.x * N1.y * V0.z - N2.x * V1.y * V0.z
                 - V0.x * V2.y * N1.z - V0.x * N2.y * V1.z - N0.x * V2.y * V1.z - V1.x * V0.y * N2.z - V1.x * N0.y * V2.z - N1.x * V0.y * V2.z;

    return 3.*t*t * a0 + 2.*t * a1 + a2;
}

double App::findCubicRoot(double a0, double a1, double a2, double a3){
    double p = -a1 / (3. * a0);
    double r = a2 / (3. * a0);

    double q = p*p*p + (a1 * a2 - 3. * a0 * a3)/(6. * a0 * a0);

    double t = std::pow(q + sqrt(q*q + std::pow(r - p * p,3)),1/3.) +  std::pow(q - sqrt(q * q + std::pow(r - p * p,3)),1/3.) + p;

    return t;
}

Vector3 App::findTriangleBarycentricCoordinates(Vector3 P, Vector3 T0, Vector3 T1, Vector3 T2){
    Vector3 u = T1 - T0;
    Vector3 v = T2 - T0;

    Vector3 n = u.cross(v);

    float oneOver4ASquared = 1.0f / n.dot(n);

    Vector3 w = P - T0;

	float b2 = u.cross(w).dot(n) * oneOver4ASquared;
	float b1 = w.cross(v).dot(n) * oneOver4ASquared;

    return Vector3(1.f - b2 - b1, b1, b2);
}
