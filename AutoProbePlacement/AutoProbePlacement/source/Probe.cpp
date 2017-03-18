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


void Probe::initProbeCoefficients(TProbeCoefficients& coeffs)
{
	int NumCoeffs = 9;
	for (int i = 0; i < NumCoeffs; ++i)
	{
		coeffs.push_back(G3D::Vector3(0, 0, 0));
	}
}
CoeffGradients CreateCoeffGradients(int NumCoefficients)
{
	CoeffGradients coeffGrad;
	for (int coeff = 0; coeff < NumCoefficients; ++coeff)
	{
		coeffGrad.append(Array<Vector3>());

	}
	for (int coeff = 0; coeff < NumCoefficients; ++coeff)
	{
		for (int channel = 0; channel < 3; ++channel)
		{
			coeffGrad[coeff].append(Vector3());
		}
	}
	return coeffGrad;
}


Probe::Probe()
{
	
}


Probe::Probe(int index, String probeStructurePath)
{
	this->probeStructurePath = probeStructurePath;
	this->index = index;
	texture = NULL;
	irr_texture = NULL;

	initProbeCoefficients(this->coeffs);
	this->coeffGradients = CreateCoeffGradients(9);
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
		case EResource::Probes_X_NEG:
			return (this->probeStructurePath) + "/Probes/Probe_" + String(std::to_string(index).c_str()) + "_0.png";
		case EResource::Probes_X_POS:
			return (this->probeStructurePath) + "/Probes/Probe_" + String(std::to_string(index).c_str()) + "_1.png";
		case EResource::Probes_Y_NEG:
			return (this->probeStructurePath) + "/Probes/Probe_" + String(std::to_string(index).c_str()) + "_2.png";
		case EResource::Probes_Y_POS:
			return (this->probeStructurePath) + "/Probes/Probe_" + String(std::to_string(index).c_str()) + "_3.png";
		case EResource::Probes_Z_NEG:
			return (this->probeStructurePath) + "/Probes/Probe_" + String(std::to_string(index).c_str()) + "_4.png";
		case EResource::Probes_Z_POS:
			return (this->probeStructurePath) + "/Probes/Probe_" + String(std::to_string(index).c_str()) + "_5.png";
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

void Probe::setPosition(G3D::Vector3& pos)
{
    this->frame.translation = pos;
}



void Probe::computeCoefficients(std::shared_ptr<G3D::Image> probeTexture, TProbeCoefficients& coeffs)
{
	int NumCoeffs = 9;
	int probeTextureWidth = probeTexture->width();
	int probeTextureHeight = probeTexture->height();

	for (int h = 0; h < probeTextureHeight; ++h)
	{
		for (int w = 0; w < probeTextureWidth; ++w)
		{
			Vector2 UV = Transform::ijToUV(Vector2((float)w, (float)h), probeTextureWidth, probeTextureHeight);

			Vector2 PT = Transform::UVtoPT(UV);

			float domega = (PI / (float)probeTextureHeight) * (2.f * PI / (float)probeTextureWidth) * sin(PT.y);

			Color3 radiance;
			probeTexture->get(Point2int32(w, h), radiance);

			if (radiance == Color3(0, 0, 0))
			{
				continue;
			}

			for (int k = 0; k < NumCoeffs; ++k)
			{
				// Compute contribution to the gradient coefficients
				std::pair<int, int> lm = SH::kToLM(k);
				Vector3 texelWSVector = Transform::sphericalToCartesian(PT, 1.0f);


				//float shFunctionEvaluation = SH::SHxyz_yup(lm.first, lm.second, directionVector);

				//Vector3 shFunctionEvaluationGrad = Vector3(shGradient[AXIS_X][k], shGradient[AXIS_Y][k], shGradient[AXIS_Z][k]);

				//Vector3 integralFactor0Term0 = shFunctionEvaluationGrad * g;
				//Vector3 integralFactor0Term1 = gGrad * shFunctionEvaluation;

				//Vector3 integralFactor0 = integralFactor0Term0 + integralFactor0Term1;

				//Array<Vector3> integralTerm;
				//for (int channel = 0; channel < 3; ++channel)
				//{
				//	integralTerm.append(integralFactor0 * radiance[channel] * dA);
				//}

				//for (int channel = 0; channel < 3; ++channel)
				//{
				//	tempCoeffGradients[k][channel] = tempCoeffGradients[k][channel] + integralTerm[channel];
				//}

				////if ((k == 0) && (h == 15) && (w == 0))
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
				coeffs[k] += toAdd;
			}
		}
	}
}


TProbeCoefficients subtractAndDivide(TProbeCoefficients& c0, TProbeCoefficients& c1, float divider)
{
	TProbeCoefficients diff;
	Probe::initProbeCoefficients(diff);

	for (int i = 0; i < diff.size(); ++i)
	{
		diff[i] = (c1[i] - c0[i]) / divider;
	}

	return diff;
}

CoeffGradients Probe::computeProbeCoeffGradients(float gradientDisplacement)
{
	int NumCoeffs = 9;
	CoeffGradients toReturn = CreateCoeffGradients(NumCoeffs);

	String paths[6];
	paths[0] = buildPath(EResource::Probes_X_NEG);
	paths[1] = buildPath(EResource::Probes_X_POS);
	paths[2] = buildPath(EResource::Probes_Y_NEG);
	paths[3] = buildPath(EResource::Probes_Y_POS);
	paths[4] = buildPath(EResource::Probes_Z_NEG);
	paths[5] = buildPath(EResource::Probes_Z_POS);

	try
	{
		std::shared_ptr<G3D::Image> textures[6];
		textures[0] = Image::fromFile(paths[0]);
		textures[1] = Image::fromFile(paths[1]);
		textures[2] = Image::fromFile(paths[2]);
		textures[3] = Image::fromFile(paths[3]);
		textures[4] = Image::fromFile(paths[4]);
		textures[5] = Image::fromFile(paths[5]);

		TProbeCoefficients dCoeffs[6];
		for (int c = 0; c < 6; ++c)
		{
			initProbeCoefficients(dCoeffs[c]);
			computeCoefficients(textures[c], dCoeffs[c]);
		}

		TProbeCoefficients gradients[3];

		float distance = gradientDisplacement;
		gradients[0] = subtractAndDivide(dCoeffs[0], dCoeffs[1], distance);
		gradients[1] = subtractAndDivide(dCoeffs[2], dCoeffs[3], distance);
		gradients[2] = subtractAndDivide(dCoeffs[4], dCoeffs[5], distance);

		for (int band = 0; band < NumCoeffs; ++band)
		{
			for (int channel = 0; channel < 3; ++channel)
			{
				for (int axis = 0; axis < 3; ++axis)
				{
					toReturn[band][channel][axis] = gradients[axis][band][channel];
				}
			}
		}
	}
	catch (...)
	{

	}
	

	return toReturn;
}

void Probe::computeCoefficientsFromTexture(bool alsoSet, bool computeGradients, float gradientDisplacement)
{
	int NumCoeffs = 9;

	TProbeCoefficients tempCoeffs;
	CoeffGradients tempCoeffGradients = CreateCoeffGradients(NumCoeffs);

	initProbeCoefficients(tempCoeffs);

	String texturePath = buildPath(EResource::Probes);
	String positionPath = buildPath(EResource::Positions);
	String normalPath = buildPath(EResource::Normals);

	try
	{
		std::shared_ptr<G3D::Image> probeTexture = Image::fromFile(texturePath);
		computeCoefficients(probeTexture, tempCoeffs);
	}
	catch (G3D::Image::Error e)
	{
		debugPrintf("Exception when loading probe texture (filename='%s'): %s\n", e.filename, e.reason);
	}

	if (alsoSet)
	{
		this->coeffs = tempCoeffs;
		if (computeGradients)
		{
			this->coeffGradients = computeProbeCoeffGradients(gradientDisplacement);
		}
	}
}

void Probe::saveCoefficients()
{
	String coeffFilePath = buildPath(EResource::Coefficients);
	String gradientFilePath = buildPath(EResource::CoefficientGradients);

	std::fstream coeffFile;
	std::fstream gradientFile;
	coeffFile.open(coeffFilePath.c_str(), std::ios::out);
	gradientFile.open(gradientFilePath.c_str(), std::ios::out);
	coeffFile.precision(20);
	gradientFile.precision(20);

	char channels[3] = { 'R', 'G', 'B' };

	for (int i = 0; i < coeffs.size(); ++i)
	{
		gradientFile << "// Coeff #" << i << std::endl;

		coeffFile << coeffs[i][0] << std::endl;
		coeffFile << coeffs[i][1] << std::endl;
		coeffFile << coeffs[i][2] << std::endl;

		for (int channel = 0; channel < 3; ++channel)
		{
			gradientFile << "// Channel " << channels[channel] << std::endl;
			for (int axis = AXIS_X; axis <= AXIS_Z; ++axis)
			{
				gradientFile << coeffGradients[i][channel][axis] << std::endl;
			}
		}
	}

	coeffFile.close();
	gradientFile.close();
}

void Probe::reconstructSH(const G3D::Vector3& normal)
{
	for (int i = 0; i < coeffs.size(); ++i)
	{

	}
}