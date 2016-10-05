#include <Windows.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include "Helpers.h"

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

bool runPythonScriptFromDataFiles(std::string scriptName, std::string args)
{
	std::stringstream ss;

	ss << "cmd /c \"cd C:\\git\\AutoProbePlacement\\AutoProbePlacement\\data-files\\scripts && C:\\Users\\Joel\\Anaconda2\\python.exe ";
	//ss << "cmd /c \"cd C:\\git\\AutoProbePlacement\\AutoProbePlacement\\data-files\\scripts && python ";
	ss << scriptName << " ";
	ss << args << "\"";

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