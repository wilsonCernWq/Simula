#pragma once

#include <iostream>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"
#include "xmlNodeDefinition_Action.h"

namespace simula {

	///------------------------------------------------------
	/// Occurrence
	struct Occurrence_Core {
	public:
		UINT gIdx; ///< global index (used in KMC lookup)
		UINT rIdx; ///< relative (local) index
		UINT iid_action, ///< initial action index
			 fid_action; ///< final action index
	public:
		UINT occur_min, occur_max;
	};
	class Occurrence : public NodeImplementation<Occurrence_Core> {
	protected:
		std::vector<Action> actions;
		CHAR component_id_char;
	public:
		void attribute(CHAR str, CHAR content) {
			if (str.compare("min") == 0) {
				core->occur_min = stoi(content);
			}
			else if (str.compare("max") == 0) {
				core->occur_max = stoi(content);
			}
			else if (str.compare("cid") == 0) {
				component_id_char = content;
			}
			else {
				std::cerr << "Error: undefined attribute found " << str << std::endl;
				throw 0;
			}
		}
		Node* value(CHAR str) {
			if (str.compare("action") == 0) {
				actions.emplace_back();
				return &actions.back();
			}
			else {
				std::cerr << "Error: undefined node found " << str << "n";
				throw 0; return nullptr;
			}
		}
		// TODO add function to deep copy data into core type
		// TODO add update function to link the array counter for action
	};

}