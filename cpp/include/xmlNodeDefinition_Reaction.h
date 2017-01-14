#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"
#include "xmlNodeDefinition_Spot.h"
#include "xmlNodeDefinition_Occurrence.h"

namespace simula {
	///------------------------------------------------------
	/// Reaction
	struct Reaction_Core {
		UINT gIdx, ///< global index (used in KMC lookup)
			rIdx; ///< relative (local) index
		UINT iid_spot,
			fid_spot;
		UINT iid_occurrence,
			fid_occurrence;
		/// --- fields from input
		REAL energy; ///< reaction energ
		INT state_i; ///< initial state
		INT state_f; ///< final state
	};
	std::vector<Reaction_Core> vector_reactions;
	class Reaction : public Node {
		std::shared_ptr<Reaction_Core> core;
		std::vector<Spot> spots;
		std::vector<Occurrence> occurrences;
	public:
		Reaction() {
			vector_reactions.emplace_back();
			core = std::make_shared<Reaction_Core>(vector_reactions.back());
		}
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
		std::shared_ptr<Node> value(CHAR str) {
			if (str.compare("spot") == 0) {
				auto node = std::make_shared<Spot>();
				spots.push_back(*node);
				return node;
			}
			else if (str.compare("occurrence") == 0) {
				auto node = std::make_shared<Occurrence>();
				occurrences.push_back(*node);
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