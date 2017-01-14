#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"

namespace simula {
	///------------------------------------------------------
	///< Action
	struct Action_Core {
		UINT gId; ///< global index (used in KMC lookup)
		UINT rId; ///< relative (local) index
				  /// --- fields from input
		UINT component_id;
		INT state_i; ///< initial state
		INT state_f; ///< final state
	};
	std::vector<Action_Core> vector_actions;
	class Action : public Node {
		std::shared_ptr<Action_Core> core;
		CHAR component_id_char;
	public:
		Action() {
			vector_actions.emplace_back();
			core = std::make_shared<Action_Core>(vector_actions.back());
		}
		void attribute(CHAR str, CHAR content) {
			if (str.compare("cid") == 0) {
				component_id_char = content;
			}
			else if (str.compare("initial") == 0) {
				core->state_i = stoi(content);
			}
			else if (str.compare("final") == 0) {
				core->state_f = stoi(content);
			}
			else {
				std::cerr << "Error: undefined attribute " << str << std::endl;
				throw 0;
			}
		}
		/// TODO add function to deep copy data into core type
		/// TODO function to update component id into index
	};
}