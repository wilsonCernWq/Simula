#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"
#include "xmlNodeDefinition_Position.h"

namespace simula {
	///------------------------------------------------------
	/// Molecule
	struct Molecule_Core {
		UINT gid;
		UINT size;
		UINT iid_position,
			fid_position;
		/// --- fields from input
		CHAR id_char;
		INT x, y, z, d;
	};
	std::vector<Molecule_Core> vector_molecules;
	class Molecule : public Node {
		std::shared_ptr<Molecule_Core> core;
		std::vector<Position> positions;
	public:
		Molecule() {
			vector_molecules.emplace_back();
			core = std::make_shared<Molecule_Core>(vector_molecules.back());
		}
		std::shared_ptr<Node> value(std::string str) {
			if (str.compare("position") == 0) {
				auto node = std::make_shared<Position>();
				this->positions.push_back(*node);
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