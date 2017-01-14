#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <../../utils/rapidxml_1_13/rapidxml/rapidxml.hpp>

#include "simulaBaseType.h"
#include "xmlNodeDefinition_Document.h"

namespace simula {

	/** Function for reading XML input files
	*/
	void readnode(rapidxml::xml_node<>* n, INT level, Node* parent)
	{
		using namespace rapidxml;
		using namespace std;
		typedef xml_attribute<> xmlattr;
		typedef xml_node<> xmlnode;
		std::string indent = "    ";

		/// TODO define objects here
		Node* child = parent->value(string(n->name()));

		///< node information
		for (INT i = 0; i < level; ++i) { cout << indent; } ///< level indentation
		cout << " << " << n->name() << " >> ";
		cout << " " << n->value() << "\n";

		///< read attributes
		for (xmlattr* a = n->first_attribute(); a; a = a->next_attribute()) {
			for (INT i = 0; i < level; ++i) { cout << indent; } ///< level indentation
																///< format: [...] attribute: xxx value: xxx
			cout << " ** " << a->name();
			cout << " = " << a->value() << "\n";

			child->attribute(a->name(), a->value());
		}


		///< read children
		for (xmlnode* c = n->first_node(); c; c = c->next_sibling()) {
			readnode(c, level + 1, child);
		}

	}

	void parseInput(std::string& text) {

		/// read xml using rapidxml
		rapidxml::xml_document<> doc;
		doc.parse<0>(&text[0]);

		/// print out data
		auto rNode = doc.first_node();

		/// for all sub nodes
		auto root = new simula::Document();
		for (auto c = rNode->first_node(); c; c = c->next_sibling()) {
			simula::readnode(c, 0, root);
		}
		std::cout << "done here \n";
		simula::Action::update();
	}

}
