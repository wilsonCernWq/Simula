#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"
#include "xmlNodeDefinition_Action.h"

namespace simula {
	///------------------------------------------------------
	/// Occurrence
	struct Occurrence_Core {
		UINT gIdx; ///< global index (used in KMC lookup)
		UINT rIdx; ///< relative (local) index
		UINT iid_action, ///< initial action index
			fid_action; ///< final action index
						/// --- fields from input
		UINT occur_min, occur_max;
	};
	std::vector<Occurrence_Core> vector_occurrences;
	class Occurrence : public Node {
		std::shared_ptr<Occurrence_Core> core;
		std::vector<Action> actions;
	public:
		Occurrence() {
			vector_occurrences.emplace_back();
			core = std::make_shared<Occurrence_Core>(vector_occurrences.back());
		}
		void attribute(CHAR str, CHAR content) {
			if (str.compare("min") == 0) {
				core->occur_min = stoi(content);
			}
			else if (str.compare("max") == 0) {
				core->occur_max = stoi(content);
			}
			else {
				std::cerr << "Error: undefined attribute found " << str << std::endl;
				throw 0;
			}
		}
		std::shared_ptr<Node> value(CHAR str) {
			if (str.compare("action") == 0) {
				auto node = std::make_shared<Action>();
				actions.push_back(*node);
				return node;
			}
			else {
				std::cerr << "Error: undefined node found " << str << "n";
				throw 0;
				return nullptr;
			}
		}
		// TODO add function to deep copy data into core type
		// TODO add update function to link the array counter for action
	};
}