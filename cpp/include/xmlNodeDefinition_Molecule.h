#pragma once

#include <iostream>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"
#include "xmlNodeDefinition_Position.h"

namespace simula {

	///------------------------------------------------------
	/// Molecule
	struct Molecule_Core {
	public:
		UINT gid;
		UINT size;
		UINT iid_position, fid_position;
		INT x, y, z, d;
	public:
		CHAR id_char;
	};
	class Molecule : public NodeImplementation<Molecule_Core> {
	protected:
		std::vector<Position> positions;
	public:
		void attribute(std::string str, std::string content) {
			if (str.compare("id") == 0) {
				core->id_char = content;
			}
			else {
				std::cerr << "Error: undefined attribute " << str << std::endl;
				throw 0;
			}
		}
		Node* value(std::string str) {
			if (str.compare("position") == 0) {
				positions.emplace_back();
				return &positions.back();
			}
			else {
				std::cerr << "Error: undefined node found " << str << "n";
				throw 0; return nullptr;
			}
		}
	};

}