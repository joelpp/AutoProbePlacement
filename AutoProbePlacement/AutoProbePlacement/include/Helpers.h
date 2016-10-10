#pragma once
#include <string>
#include "G3D\G3DAll.h"
#define PI 3.141592654f

bool runCommand(std::string command);

bool runPythonScriptFromDataFiles(std::string scriptName, std::string args, bool showOutput);

bool createFolder(const char* name);

void createEmptyFile(const char* name);

double NDotOmegaCoeff(int l);

float phongCoeffs(int l, float r);

Color3 sRGBtoRGB(const Color3& source);
//G3D::Array<G3D::String> getFoldersInFolder(const G3D::String& path);
