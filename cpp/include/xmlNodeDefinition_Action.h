#pragma once

#include <iostream>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"

namespace simula {

	///------------------------------------------------------
	///< Action
	struct Action_Core {
	public:
		UINT gId; ///< global index (used in KMC lookup)
	public:
		UINT component_id;
		INT state_i; ///< initial state
		INT state_f; ///< final state
	};
	class Action : public NodeImplementation<Action_Core> {
	protected:
		CHAR component_id_char;
	public:
		virtual void attribute(CHAR str, CHAR content) {
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
		/// update the information inside 
		void static update() {
			for (auto i = 0; i != Action::list.size(); ++i) {
				auto it = Action::list[i];
				Action::list[i].gId = i;
				std::cout << "# " << i << " " << it.state_i << " " << it.state_f << std::endl;
			}
		}
	};

}