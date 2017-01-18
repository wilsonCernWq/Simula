#pragma once

#include <iostream>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"
#include "xmlNodeDefinition_Reaction.h"

namespace simula {

	///------------------------------------------------------
	/// Component
	struct Component_Core {
	public:
		UINT gId;
		UINT iid_reaction, fid_reaction;
	public:
		UINT component_id;      ///< component index
		CHAR component_id_char; ///< component index defined
		INT  bgcolor;           ///< initial color (backgrounf color)
	};
	class Component : public NodeImplementation<Component_Core> {
	protected:
		std::vector<Reaction> reactions;
	public:
		void attribute(std::string str, std::string content) {
			if (str.compare("cid") == 0) {
				core->component_id_char = content;
			}
			else if (str.compare("bgcolor") == 0) {
				core->bgcolor = stoi(content);
			}
			else {
				std::cerr << "Error: undefined attribute " << str << std::endl;
				throw 0;
			}
		}
		Node* value(CHAR str) {
			if (str.compare("reaction") == 0) {
				reactions.emplace_back();
				return &reactions.back();
			}
			else {
				std::cerr << "Error: undefined node found " << str << "n";
				throw 0; return nullptr;
			}
		}
		// update list information after definition
		void update() {
			for (auto i = 0; i < Component::list.size(); ++i) {
				auto it = Component::list[i];
				Component::list[i].gId = i;

			}
		}
	};

}
