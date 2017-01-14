#pragma once

#include <iostream>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"

namespace simula {

	///------------------------------------------------------
	/// Position
	struct Position_Core {
		UINT id_component;
		INT  rx, ry, rz, rd;
	};
	class Position : public NodeImplementation<Position_Core> {
	protected:
		CHAR component_id_char;
	public:
		void attribute(std::string str, std::string content) {
			if (str.compare("component") == 0) {
				component_id_char = content;
			}
			else if (str.compare("rx") == 0) {
				core->rx = stoi(content);
			}
			else if (str.compare("ry") == 0) {
				core->ry = stoi(content);
			}
			else if (str.compare("rz") == 0) {
				core->rz = stoi(content);
			}
			else if (str.compare("rd") == 0) {
				core->rd = stoi(content);
			}
			else {
				std::cerr << "Error: undefined attribute " << str << std::endl;
				throw 0;
			}
		}
	};
}
