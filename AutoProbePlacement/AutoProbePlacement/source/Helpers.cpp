#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include "Helpers.h"
#include "SH.h"

String generateFolderNameBaseAnySuffix(const String& prefix) 
{
    Array<String> exist;
    String templat = prefix;
    FileSystem::getDirectories(templat + "*", exist, true);

    int num = 0;
    String result;
    templat += "%d";
    bool done;
    do {
        result = format(templat.c_str(), num);
        ++num;
        done = true;
        for (int f = 0; f < exist.length(); ++f) {
            if (beginsWith(exist[f], result)) {
                done = false;
                break;
            }
        }
    } while (!done);
    return result;
}

int folderCount(const String& path)
{
    Array<String> exist;

    FileSystem::getDirectories(path + "*", exist, true);
    return exist.size();
}

bool runCommand(std::string command)
{
	//std::stringstream args;
	//args << "cmd /c \"" << command <<" \"";

	STARTUPINFOA si;
	PROCESS_INFORMATION pi;

	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));
	std::string comm = "C:\\Windows\\winsxs\\amd64_microsoft-windows-commandprompt_31bf3856ad364e35_6.1.7601.17514_none_e932cc2c30fc13b0\\cmd.exe";
	LPCSTR sw = comm.c_str();
	LPSTR arg2 = const_cast<char *>(command.c_str());
	bool result = CreateProcessA(sw, arg2, NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi);

	if (result == false)
	{
		DWORD err = GetLastError();
		std::cout << "ERROR:" << err << std::endl;
	}
	else if (result == true)
	{
		WaitForSingleObject(pi.hProcess, INFINITE);
		CloseHandle(pi.hProcess);
		CloseHandle(pi.hThread);
	}

	return result;
}

bool runPythonScriptFromDataFiles(std::string scriptName, std::string args, bool showOutput)
{
	std::stringstream ss;

	//ss << "cmd /c \"cd C:\\git\\AutoProbePlacement\\AutoProbePlacement\\data-files\\scripts && C:\\Users\\Joel\\Anaconda2\\python.exe ";
	ss << "cmd /c \"cd C:\\git\\AutoProbePlacement\\AutoProbePlacement\\data-files\\scripts && python ";
	ss << scriptName << " ";
	ss << args << "\"";

	if (!showOutput)
	{
		ss << " > nul";
	}

	std::cout << ss.str();
	return runCommand(ss.str());
}

bool createFolder(const char* name)
{
	CreateDirectoryA(name, NULL);
	return ERROR_ALREADY_EXISTS == GetLastError();
}

void createEmptyFile(const char* name)
{
	std::fstream file(name, std::fstream::out);
}
//
//
//G3D::Array<G3D::String> getFoldersInFolder(const G3D::String& path)
//{
//	G3D::Array<G3D::String> folderList;
//	FileSystem::ListSettings ls;
//	ls.includeParentPath = false;
//	ls.recursive = false;
//
//	FileSystem::list(path + "/*", folderList, ls);
//
//	return folderList;
//}


float phongCoeffs(int l, float r)
{
    if (l == 0)
    {
        return PI;
    }
    else if (l == 1)
    {
        return PI * (1.0f + r) / (2.0f + r);
    }
    else if (l == 2)
    {
        return PI * r / (3.0f + r);
    }
    return 0;
}

double NDotOmegaCoeff(int l)
{
    if (l == 1)
    {
        return 2 * PI / 3;
    }
    else if (l % 2 == 0)
    {
        int a = (int)pow(-1, l / 2.f - 1);
        int b = SH::factorial(l);
        return (double)((2 * PI * a * b) / ((l + 2) * (l + 1) * ((int)pow(2, l)) * SH::factorial((int)pow(((float)l) / 2.f, 2))));
    }
    else
    {
        return 0;
    }
}

//def sRGBToRGBVal(val) :
//    if (val < 0.04045) :
//        return val / 12.92;
//    else :
//        return ((val + 0.055) / 1.055)**2.4;
//
//def sRGBToRGB(srgbColor) :
//    rgb = [];
//for x in srgbColor :
//rgb.push(sRGBToRGBVal(x))
//
//return rgb;

float sRGBtoRGB(float source)
{
    float c = 0;

    if (source < 0.04045)
    {
        return source / 12.92;
    }
    else
    {
        return pow(((source + 0.055) / 1.055), 2.4);
    }
}
Color3 sRGBtoRGB(const Color3& source)
{
    return Color3(sRGBtoRGB(source.r), sRGBtoRGB(source.g), sRGBtoRGB(source.b));
}