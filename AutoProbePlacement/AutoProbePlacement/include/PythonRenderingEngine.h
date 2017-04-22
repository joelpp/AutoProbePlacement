#pragma once

#include <string>


class PythonRenderingEngine
{
public:
	PythonRenderingEngine();

	PythonRenderingEngine(bool running);

	void setControlFilePath();

	void start(const char* sceneName);

	void stop();

	bool isRunning() { return m_running; }
private:
	bool m_running;
	std::string m_ControlFilePath;


};