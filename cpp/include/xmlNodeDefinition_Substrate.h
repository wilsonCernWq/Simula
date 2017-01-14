#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"

namespace simula {
	/** Definition of a substrate
	* Matrices are stored in row-major order:
	* M(row, col) = *(M.elements + row * M.width + col)
	*/
	class Substrate : public Node {
		INT width;
		INT height;
		float* elements;
	public:
		void attribute(std::string str, std::string content) {
			if (str.compare("width") == 0) {
				this->width = stoi(content);
			}
			else if (str.compare("height") == 0) {
				this->height = stoi(content);
			}
			else {
				std::cerr << "undefined attribute " << str << std::endl;
				throw 0;
			}
		}
	};
}
