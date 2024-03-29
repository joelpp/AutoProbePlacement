#pragma once
#include <string>
#include "G3D\G3DAll.h"
#define PI 3.141592654f
#define DEAV3(x) debugPrintf(#x); for(int num = 0; num < x.size(); num++){ debugPrintf(", [%d]: (%s)\n", num, x[num].toString().c_str()); debugPrintf("\n");}
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

#define EMPTY_THROW(v) try { #v; } catch (std::exception e) { };

String generateFolderNameBaseAnySuffix(const String& prefix);

bool runCommand(std::string command, bool waitForCompletion);
bool runShellCommand(std::string command, bool showOutput, bool waitForCompletion);

bool runPythonScriptFromDataFiles(std::string scriptName, std::string args, bool showOutput, bool waitForCompletion);

bool createFolder(const char* name);
bool createFolder(String& name);

std::fstream createEmptyFile(const char* name);
std::fstream createEmptyFile(String name);

void copyDir(const char* srcPath, const char* dstPath);
void copyDir(const String& srcPath, const String& dstPath);

void copyFile(const char* srcPath, const char* dstPath);
void copyFile(const String& srcPath, const String& dstPath);

int folderCount(const String& path);

double NDotOmegaCoeff(int l);

float phongCoeffs(int l, float r);

Color3 sRGBtoRGB(const Color3& source);

void dumpToFile(std::fstream& file, const Array<Vector3>& arr);
void dumpToFile(std::fstream& file, const Array<float>& arr);

void dumpZerosToFile(std::fstream& file, int amt);

Vector3 StringToVector3(const G3D::String& s);
Vector3 StringToVector3(const std::string& s);

template<typename T>
String g3dString(T s)
{
    return String(std::to_string(s).c_str());
};

void popNotification(std::string title, std::string message, int timeInSeconds);

void normalize(G3D::Vector3& v);

std::vector<float> readValuesFromFlatFile(const char* fileName);

void writeValuesToFlatFile(const char* fileName, std::vector<float>& values);

FILETIME getFileLastModifiedTime(const char* fileName);

bool isLaterFileTime(FILETIME tested, FILETIME reference);