#pragma once

#include <iostream>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"
#include "xmlNodeDefinition_Spot.h"
#include "xmlNodeDefinition_Occurrence.h"

namespace simula {

	///------------------------------------------------------
	/// Reaction
	struct Reaction_Core {
	public:
		UINT gIdx; ///< global index (used in KMC lookup)
		UINT iid_spot, fid_spot;
		UINT iid_occurrence, fid_occurrence;
	public:
		REAL energy; ///< reaction energ
		INT state_i; ///< initial state
		INT state_f; ///< final state
	};
	class Reaction : public NodeImplementation<Reaction_Core> {
	protected:
		std::vector<Spot> spots;
		std::vector<Occurrence> occurrences;
	public:
		void attribute(std::string str, std::string content) {
			if (str.compare("energy") == 0) {
				core->energy = stof(content);
			}
			else if (str.compare("initial") == 0) {
				core->state_i = stoi(content);
			}
			else if (str.compare("final") == 0) {
				core->state_f = stoi(content);
			}
			else {
				std::cerr << "Error: undefined attribute found " << str << std::endl;
				throw 0;
			}
		}
		Node* value(CHAR str) {
			if (str.compare("spot") == 0) {
				spots.emplace_back();
				return &spots.back();
			}
			else if (str.compare("occurrence") == 0) {
				occurrences.emplace_back();
				return &occurrences.back();
			}
			else {
				std::cerr << "Error: undefined node found " << str << "n";
				throw 0; return nullptr;
			}
		}
	};

}