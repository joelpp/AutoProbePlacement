#include "Transform.h"
#include "Probe.h"
#include "SH.h"
#include <G3D/G3DAll.h>
#include <fstream>
#include <sstream>

#define PI 3.141592654f
#define AXIS_X 0
#define AXIS_Y 1
#define AXIS_Z 2

Probe::Probe()
{
	
}


Probe::Probe(int index, String probeStructurePath)
{
	this->probeStructurePath = probeStructurePath;
	this->index = index;
	texture = NULL;
	irr_texture = NULL;

}

String Probe::buildPath(EResource res)
{
	// TODO: clean up this mess...
	switch(res)
	{
		case EResource::Probes: 
								return (this->probeStructurePath) + "/Probes/Probe_" + String(std::to_string(index).c_str()) + ".png";
		case EResource::Normals: 
								return (this->probeStructurePath) + "/Normals/Normal_" + String(std::to_string(index).c_str()) + ".exr";
		case EResource::Positions: 
								return (this->probeStructurePath) + "/Positions/Position_" + String(std::to_string(index).c_str()) + ".exr";
		case EResource::Coefficients: 
								return (this->probeStructurePath) + "/Coefficients/Coefficients_" + String(std::to_string(index).c_str()) + ".txt";
		case EResource::CoefficientGradients: 
								return (this->probeStructurePath) + "/Coefficients/CoefficientsGradients_" + String(std::to_string(index).c_str()) + ".txt";
	};
	return "INVALID RESOURCE TYPE";
}


shared_ptr<Texture> Probe::getTexture(int texID)
{
	if (bNeedsUpdate)
	{
		return Texture::fromFile("<white>");
	}

	switch(texID)
	{
		case 0: if (!texture)
				{
					String filename = buildPath(EResource::Probes);

					try
					{
						shared_ptr<Image> img = Image::fromFile(filename);
						img->flipHorizontal();
						const G3D::ImageFormat* format = G3D::ImageFormat::SRGB8();;
						texture = Texture::fromImage(G3D::String(filename), img, format);
					}
					catch(...)
					{
						texture = Texture::fromFile("<white>");
					}
				}
				return texture;
				break;
			
		case 1: return irr_texture;
				break;
		default:
				return NULL;
				break;
	}
}


shared_ptr<ProbeManipulator> Probe::getManipulator()
{
	return manipulator;
}


Point3 Probe::getPosition()
{
	return this->frame.translation;
}


CoeffGradients CreateCoeffGradients(int NumCoefficients)
{
	CoeffGradients coeffGrad;
	for (int coeff = 0 ; coeff < NumCoefficients; ++coeff)
	{
		coeffGrad.append(Array<Vector3>());

	}
	for (int coeff = 0 ; coeff < NumCoefficients; ++coeff)
	{	
		for (int channel = 0; channel < 3; ++channel)
		{
			coeffGrad[coeff].append(Vector3());
		}
	}
	return coeffGrad;
}

void Probe::computeCoefficientsFromTexture(bool alsoSet)
{
	int NumCoeffs = 9;

	G3D::Array<G3D::Vector3> tempCoeffs;
	CoeffGradients tempCoeffGradients = CreateCoeffGradients(NumCoeffs);

	for (int i = 0; i < NumCoeffs; ++i)
	{
		tempCoeffs.push_back(G3D::Vector3(0,0,0));
	}

	String texturePath = buildPath(EResource::Probes);
	String positionPath = buildPath(EResource::Positions);
	String normalPath = buildPath(EResource::Normals);

	std::shared_ptr<G3D::Image> probeTexture = Image::fromFile(texturePath);
	std::shared_ptr<G3D::Image> positionTexture = Image::fromFile(positionPath);
	std::shared_ptr<G3D::Image> normalTexture = Image::fromFile(normalPath);

	int probeTextureWidth = probeTexture->width();
	int probeTextureHeight = probeTexture->height();

	for (int h = 0; h < probeTextureHeight; ++h)
	{
		for (int w = 0; w < probeTextureWidth; ++w)
		{
			Vector2 UV = Transform::ijToUV(Vector2( (float) w, (float) h ), probeTextureWidth, probeTextureHeight);

			Vector2 PT = Transform::UVtoPT(UV);

			Vector3 texelWSVector = Transform::sphericalToCartesian(PT, 1.0f);

			float domega = (PI / (float)probeTextureHeight) * (2.f * PI / (float)probeTextureWidth) * sin(PT.y);
			
			Color3 radiance;
			probeTexture->get(Point2int32(w, h), radiance);
            
			Color3 cPosition;
			positionTexture->get(Point2int32(w, h), cPosition);
			Vector3 vPosition = 25.f * ((Vector3(cPosition) * 2) - Vector3(1,1,1));

			Color3 cNormal;
			normalTexture->get(Point2int32(w, h), cNormal);
			Vector3 vNormal = (Vector3(cNormal) * 2) - Vector3(1,1,1);

			Vector3 s = (vPosition - this->getPosition());
			float sNorm = s.magnitude();

			Vector3 directionVector = s.fastUnit();

			float normalDotDirectionVector = vNormal.dot(directionVector);

			float dA = domega * pow(sNorm, 2) / normalDotDirectionVector;

			float g = normalDotDirectionVector / pow(sNorm, 2);

			Vector3 gGradTerm0 = vNormal / pow(sNorm, 3);
			float gGradTerm1Factor = 3 * vNormal.dot(s) / pow(sNorm, 5);
			Vector3 gGradTerm1 = s * gGradTerm1Factor;

			Vector3 gGrad = gGradTerm0 - gGradTerm1;
			SHGradient shGradient = SH::gradients(vNormal);


			if (radiance == Color3(0, 0, 0))
			{
				continue;
			}

			for (int k = 0; k < NumCoeffs; ++k)
			{
				// Compute contribution to the gradient coefficients
				std::pair<int, int> lm = SH::kToLM(k);


				float shFunctionEvaluation = SH::SHxyz_yup(lm.first, lm.second, directionVector);

				Vector3 shFunctionEvaluationGrad = Vector3(shGradient[AXIS_X][k], shGradient[AXIS_Y][k], shGradient[AXIS_Z][k]);
				
				Vector3 integralFactor0Term0 = shFunctionEvaluationGrad * g;
				Vector3 integralFactor0Term1 = gGrad * shFunctionEvaluation;
				
				Vector3 integralFactor0 = integralFactor0Term0 + integralFactor0Term1;
				
				Array<Vector3> integralTerm;
				for (int channel = 0; channel < 3; ++channel)
				{
					integralTerm.append(integralFactor0 * radiance[channel] * dA);
				}
				
				for (int channel = 0; channel < 3; ++channel)
				{
					tempCoeffGradients[k][channel] = tempCoeffGradients[k][channel] +  integralTerm[channel];
				}

				//if ((k == 0) && (h == 15) && (w == 0))
				//{
				//	std::stringstream ss;
				//	ss << "DEBUG INFO" << std::endl;
				//	ss << vPosition << std::endl;
				//	ss << vNormal << std::endl;
				//	ss << s << std::endl;
				//	ss << normalDotDirectionVector << std::endl;
				//	ss << gGrad << std::endl;
				//	ss << shGradient[0][0] << std::endl;
				//	ss << shFunctionEvaluationGrad << std::endl;
				//	ss << integralFactor0 << std::endl;
				//	ss << integralTerm[0] << std::endl;
				//	ss << "END DEBUG INFO" << std::endl;
				//	debugPrintf("%s", ss.str().c_str());
				//}
				// Compute contribution to the SH coefficients
				float shFunctionEvaluationForTexel = SH::SHxyz_yup(lm.first, lm.second, texelWSVector); 

				Vector3 toAdd = domega * Vector3(radiance) * shFunctionEvaluationForTexel;
				tempCoeffs[k] += toAdd;

			}
		}
	}

	//this->coeffs.fastClear();
	//this->coeffs = tempCoeffs;
	//this->coeffGradients = tempCoeffGradients;

	String coeffFilePath = buildPath(EResource::Coefficients);
	String gradientFilePath = buildPath(EResource::CoefficientGradients);

	std::fstream coeffFile;
	std::fstream gradientFile;
	coeffFile.open(coeffFilePath.c_str(), std::ios::out);
	gradientFile.open(gradientFilePath.c_str(), std::ios::out);
	coeffFile.precision(20);
	gradientFile.precision(20);

	for (int i = 0 ; i < NumCoeffs; ++i)
	{
		coeffFile << tempCoeffs[i][0] << std::endl;
		coeffFile << tempCoeffs[i][1] << std::endl;
		coeffFile << tempCoeffs[i][2] << std::endl;

		for (int channel = 0; channel < 3; ++channel)
		{
			for (int axis = AXIS_X; axis <= AXIS_Z; ++axis)
			{
				gradientFile << tempCoeffGradients[i][channel][axis] << std::endl;
			}
		}
	}

	coeffFile.close();
	gradientFile.close();

	if (alsoSet)
	{
		this->coeffs = tempCoeffs;
		this->coeffGradients = tempCoeffGradients;
	}
}