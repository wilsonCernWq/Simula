///////////////////////////////////////////////////////////////////////////////
//
// tools for read xml input file
//
///////////////////////////////////////////////////////////////////////////////
#ifndef _SIMULA_READER_XML_
#define _SIMULA_READER_XML_

#include <boost/lexical_cast.hpp>      // lexical casting function
#include <rapidxml/rapidxml_utils.hpp> // read file
#include <rapidxml/rapidxml.hpp>       // xml file parser
#include <iostream>
#include <vector>
#include <map>
#include "global.h"

///////////////////////////////////////////////////////////////////////////////
// project namespace
namespace simula {
	namespace reader {

		typedef rapidxml::xml_document<>  doc_t;
		typedef rapidxml::file<>          xml_t;
		typedef rapidxml::xml_node<>      node_t;
		typedef rapidxml::xml_attribute<> attr_t;

		///////////////////////////////////////////////////////////////////////////
		//
		class NodeString {
		public:
			static NodeString* root;
		private:
			simI1 m_level;
			simString m_n; // name
			simString m_v; // value
			std::map<simString, simString>                m_attr;
			std::map<simString, std::vector<NodeString*>> m_node;
		public:
			simString value() {
				return m_v;
			}

			simString attr(simString attr_name) {
				return m_attr[attr_name];
			}

			std::vector<NodeString*> node(simString node_name) {
				return m_node[node_name];
			}

			friend NodeString* build_tree(const node_t*, const simString&, simI1);
#ifndef NDEBUG
			void print() {
				for (simI1 i = 0; i < m_level; ++i) { cout << "  "; }
				cout << "<" << m_n << "> " << m_v << " { ";
				for (auto& attr : m_attr) {
					cout << attr.first << "=" << attr.second << ", ";
				}
				cout << "}\n";
				for (auto& name : m_node) {
					for (auto& node : name.second) { node->print(); }
				}
			}
#endif
		};

		NodeString* NodeString::root = NULL;

		NodeString* build_tree(const node_t* r, const simString& node_name = "", simI1 level = 0)
		{
			NodeString* m_root = new NodeString();
			// build basic values
			m_root->m_level = level;
			m_root->m_n = node_name == "" ? simString(r->name()) : node_name;
			m_root->m_v = r->value();
			// build attributions
			for (attr_t* a = r->first_attribute(); a; a = a->next_attribute()) {
				simString name = a->name();
				if (name != "") { m_root->m_attr[name] = a->value(); }
			}
			// build sub nodes
			for (node_t* n = r->first_node(); n; n = n->next_sibling())	{
				simString name = n->name();
				if (name != "") {
					m_root->m_node[name].push_back(build_tree(n, name, level + 1));
				}
			}
			return m_root;
		}

		void make_doc(const doc_t& doc) {
			NodeString::root = build_tree(doc.first_node());
#ifndef NDEBUG
			NodeString::root->print();
#endif
		}

		///////////////////////////////////////////////////////////////////////////
		// template string parser
		template<typename T>
		T parse_str(simString str) // auto type conversion
		{
			boost::trim_if(str, boost::is_any_of(" \n\t"));
			return boost::lexical_cast<T>(str);
		}
		template<>
		simString parse_str<simString>(simString str)
		{
			return boost::trim_copy_if(str, boost::is_any_of(" \n\t"));
		}
		template<>
		simI2 parse_str<simI2>(simString str)
		{
			using boost::lexical_cast;
			using boost::is_any_of;
			using boost::token_compress_on;
			simStrVec strvec;
			simI2     result;
			// split string into string vector
			split(strvec, str, is_any_of(","), token_compress_on);
			// trim each string
			trim_if(strvec[0], is_any_of(" \n\t"));
			trim_if(strvec[1], is_any_of(" \n\t"));
			// convert string into numerical data
			result.x = lexical_cast<simI1>(strvec[0]);
			result.y = lexical_cast<simI1>(strvec[1]);
			return result;
		}

		/////////////////////////////////////////////////////////////////////////////
		// get attribute/node value
		inline const simChar* attr_value(const node_t* n, const simChar* name)
		{
			return n->first_attribute(name)->value();
		}
		inline const simChar* node_value(const node_t* n, const simChar* name)
		{
			return n->first_node(name)->value();
		}

		/////////////////////////////////////////////////////////////////////////////
		// convention function to parse node attribute
		template<typename T>
		T parse_attr(const node_t* n, const simChar* name)
		{
			return parse_str<T>(attr_value(n, name));
		}

		/////////////////////////////////////////////////////////////////////////////
		// apply function to each node
		template<typename Func>
		void for_each_node(const node_t* r, const simChar* name, Func func, const simI1 max_num = 0)
		{
			if (max_num == 0)
			{
#ifndef NDEBUG
				cout << " entering for_each under node " << r->name() << endl;
#endif
				for (node_t* n = r->first_node(name); n; n = n->next_sibling(name))
				{
#ifndef NDEBUG
					cout << "     for each debug: looping in " << r->name() << endl;
#endif
					func(n);
				}
			}
			else
			{
				simI1 max_counter = 0;
				for (node_t* n = r->first_node(name); n; n = n->next_sibling(name))
				{
					// check maximum iteration
					if (max_counter++ >= max_num) { break; }
					// execute function
					func(n);
				}
			}
#ifndef NDEBUG
			cout << " leaving for_each under node " << r->name() << endl << endl;
#endif
		}

	};
};

#endif // _SIMULA_READER_XML_
