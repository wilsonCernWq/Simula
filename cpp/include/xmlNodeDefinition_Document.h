#pragma once

#include <iostream>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"
#include "xmlNodeDefinition_Component.h"
#include "xmlNodeDefinition_Substrate.h"
#include "xmlNodeDefinition_Molecule.h"

namespace simula {

	class Document : public Node {
	public:
		std::vector<Molecule>  molecules;
		std::vector<Component> components;
		Substrate* substrate = nullptr;
	public:
		Node* value(std::string str) {
			if (str.compare("component") == 0) {
				components.emplace_back();
				return &components.back();
			}
			else if (str.compare("molecule") == 0) {
				molecules.emplace_back();
				return &molecules.back();
			}
			else if (str.compare("substrate") == 0) {
				if (this->substrate == NULL) {
					substrate = new Substrate();
				}
				else {
					std::cerr << "Error: multiple substrate definitions found\n";
					throw 0;
				}
				return substrate;
			}
			else {
				std::cerr << "Error: undefined node found " << str << "n";
				throw 0; return nullptr;
			}
		}
	};

}