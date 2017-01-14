#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"
#include "xmlNodeDefinition_Component.h"
#include "xmlNodeDefinition_Substrate.h"
#include "xmlNodeDefinition_Molecule.h"

namespace simula {
	class Document : public Node {
		std::vector<Molecule> molecules;
		std::vector<Component> components;
		std::shared_ptr<Substrate> substrate = nullptr;
		std::allocator<int>      substrate_alloc; // the default allocator for int
		std::default_delete<int> substrate_del;   // the default deleter for int
	public:
		std::shared_ptr<Node> value(std::string str) {
			if (str.compare("component") == 0) {
				auto node = std::make_shared<Component>();
				this->components.push_back(*node);
				return node;
			}
			else if (str.compare("molecule") == 0) {
				auto node = std::make_shared<Molecule>();
				this->molecules.push_back(*node);
				return node;
			}
			else if (str.compare("substrate") == 0) {
				if (this->substrate == NULL) {
					substrate = std::allocate_shared<Substrate>(substrate_alloc);
				}
				else {
					std::cerr << "Error: multiple substrate definitions found\n";
					throw 0;
				}
				return substrate;
			}
			else {
				std::cerr << "Error: undefined node found " << str << "n";
				throw 0;
			}
		}

	};
}