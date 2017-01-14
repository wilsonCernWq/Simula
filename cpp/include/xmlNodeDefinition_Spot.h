#pragma once

#include <iostream>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"

namespace simula {

	///------------------------------------------------------
	/// Spot
	struct Spot_Core {
	public:
		UINT gIdx; ///< global index (used in KMC lookup)
	public:
		UINT component_id; ///< full name for cid (converted from string to number)
		INT  rx, ry, rz, rd;
	};
	class Spot : public NodeImplementation<Spot_Core> {
		CHAR component_id_char;
	public:
		void attribute(std::string str, std::string content) {
			if (str.compare("cid") == 0) {
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
				std::cerr << "Error: undefined attribute found " << str << std::endl;
				throw 0;
			}
		}
		/// TODO function to update component id into index
	};

}