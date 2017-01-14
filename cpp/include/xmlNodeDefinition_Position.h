#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"

namespace simula {
	///------------------------------------------------------
	/// Position
	struct Position_Core {
		UINT id_component;
		UINT rx, ry, rz, rd;
	};
	std::vector<Position_Core> vector_positions;
	class Position : public Node {
		std::shared_ptr<Position_Core> core;
		CHAR component_id_char;
	public:
		Position() {
			vector_positions.emplace_back();
			core = std::make_shared<Position_Core>(vector_positions.back());
		}
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
