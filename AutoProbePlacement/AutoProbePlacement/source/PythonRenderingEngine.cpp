
#include "PythonRenderingEngine.h"
#include "Helpers.h"
#include <sstream>

PythonRenderingEngine::PythonRenderingEngine()
	: m_running(false)
{
	setControlFilePath();
}

PythonRenderingEngine::PythonRenderingEngine(bool running)
	: m_running(running)
{
	setControlFilePath();
}

void PythonRenderingEngine::setControlFilePath()
{
	if (m_ControlFilePath.empty())
	{
		m_ControlFilePath = "../data-files/scripts/optimizationSettings.txt";
	}
}

void PythonRenderingEngine::start(const char* sceneName)
{
	std::stringstream ss;
	ss << sceneName << " dummy";
	runPythonScriptFromDataFiles("onecamera_continuous.py", ss.str(), true, false);
}

void PythonRenderingEngine::stop()
{

}
