#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"

namespace simula {
	///------------------------------------------------------
	/// Spot
	struct Spot_Core {
		UINT gIdx, ///< global index (used in KMC lookup)
			rIdx; ///< relative (local) index
				  /// --- fields from input
		UINT component_id; ///< full name for cid (converted from string to number)
		INT  rx, ry, rz, rd;
	};
	std::vector<Spot_Core> vector_spots;
	class Spot : public Node {
		std::shared_ptr<Spot_Core> core;
		CHAR component_id_char;
	public:
		Spot() {
			vector_spots.emplace_back();
			core = std::make_shared<Spot_Core>(vector_spots.back());
		}
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