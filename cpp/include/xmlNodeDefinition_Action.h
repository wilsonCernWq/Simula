#pragma once

#include <iostream>
#include <vector>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Node.h"
#include "xmlNodeDefinition_Component.h"

namespace simula {

	/** 
	 * \brief Core structure definition for Action class
	 * This class defines the behaviors for all the components of molecules that surrounds
	 * the center molecule. '<cid> component_id' represents the type of the component found.
	 * '<initial> state_i' represents the initial state requirement for executing the reaction.
	 * '<final> state_f' represents the final state the component should be after the reaction.
	 */
	struct Action_Core {
	public: 
		/** fields that should be maintained by the program */
		UINT gId; ///< global index (used in KMC lookup)
	public:
		/** fields that should be read from input */
		UINT component_id;
		INT state_i; ///< initial state
		INT state_f; ///< final state
	};
	class Action : public NodeImplementation<Action_Core> {
	protected:
		///< this is the temporary storage of component index, should be converted to number later
		CHAR component_id_char; ///< since in XML file cid is a string
	public:
		Action() {
			/* call parent constructor internally before this constructor */
			core->gId = Action::list.size() - 1;
		}
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
				throw UNDEFINED_ATTRIBUTE;
			}
		}
	public:
		/// TODO add function to deep copy data into core type
		/// TODO function to update component id into index
		static UINT update_index(UINT idx) {
			if (flagDefinitionCompleted) {

			}
			else {
				std::cerr << "Error: trying to update the action list before definition is completed" << std::endl;
				throw DEFINITION_IMCOMPLETE;
			}

		}
		/// update the information inside 
		static void update() {
			for (auto it = Action::list.begin(); it != Action::list.end(); ++it) {
				std::cout << "# " << it->gId << " " << it->state_i << " " << it->state_f << std::endl;
			}
		}
	};

}