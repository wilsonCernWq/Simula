#pragma once

#include <string>

#define UNDEFINED_ATTRIBUTE   10000
#define UNDEFINED_NODE        10001
#define DEFINITION_IMCOMPLETE 10002

namespace simula {
	/// define internal types inside the program
	typedef std::string  CHAR;
	typedef int          INT;
	typedef unsigned int UINT;
	typedef double       REAL;
	typedef bool         BOOL;
	/// global flags
	BOOL flagDefinitionCompleted = false;
};
