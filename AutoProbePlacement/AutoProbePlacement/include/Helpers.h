#pragma once
#include <string>

bool runCommand(std::string command);

bool runPythonScriptFromDataFiles(std::string scriptName, std::string args);

bool createFolder(const char* name);

void createEmptyFile(const char* name);
