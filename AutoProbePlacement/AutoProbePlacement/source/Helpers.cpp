#include <Windows.h>

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
	std::string comm = "C:\\Windows\\winsxs\\amd64_microsoft-windows-commandprompt_31bf3856ad364e35_6.1.7600.16385_none_e701b864340d9016\\cmd.exe";
	LPCSTR sw = comm.c_str();
	LPSTR arg2 = const_cast<char *>(command.c_str());
	BOOL result = CreateProcessA(sw, arg2, NULL, NULL, false, CREATE_NO_WINDOW, NULL, NULL, &si, &pi);

	if (result == FALSE)
	{
		DWORD err = GetLastError();
		//debugPrintf("ERROR: %lu \n", err);
	}
	else if (result == TRUE)
	{
		WaitForSingleObject(pi.hProcess, INFINITE);
		CloseHandle(pi.hProcess);
		CloseHandle(pi.hThread);
	}

	return result;
}