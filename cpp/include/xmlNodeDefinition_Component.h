#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"
#include "xmlNodeDefinition_Reaction.h"

namespace simula {
	///------------------------------------------------------
	/// Component
	struct Component_Core {
		UINT gIdx;
		UINT iid_reaction,
			fid_reaction;
		/// --- fields from input
		UINT component_id; ///< component index
		CHAR component_id_char;
		INT  bgcolor;      ///< initial color (backgrounf color)
	};
	std::vector<Component_Core> vector_components;
	class Component : public Node {
		std::shared_ptr<Component_Core> core;
		std::vector<Reaction> reactions;
	public:
		Component() {
			vector_components.emplace_back();
			core = std::make_shared<Component_Core>(vector_components.back());
		}
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
		std::shared_ptr<Node> value(CHAR str) {
			if (str.compare("reaction") == 0) {
				auto node = std::make_shared<Reaction>();
				reactions.push_back(*node);
				return node;
			}
			else {
				std::cerr << "Error: undefined node found " << str << "n";
				throw 0;
				return nullptr;
			}
		}
	};
}
