///////////////////////////////////////////////////////////////////////////////
//
// initialization function definitions
//
///////////////////////////////////////////////////////////////////////////////
#include "init.h"
#include "substrate.h" // Substrate, Molecule
#include "reader.h"    // reader functions

///////////////////////////////////////////////////////////////////////////////
// used namespace
using namespace simula;

///////////////////////////////////////////////////////////////////////////////
// private variable/function namespace
namespace {
	/////////////////////////////////////////////////////////////////////////////
	// initialize molecules
	void init_molecule(const MoleculeType& t)
	{
		for (simI1 i = 0; i < t.eva_num(); ++i) { molecules.new_molecule(t); }
	}

	/////////////////////////////////////////////////////////////////////////////
	// initialize all molecule types
	void init_molecule_type(const reader::node_t* n)
	{
		using namespace reader;
#ifndef NDEBUG
		cout << " >>> new molecule <<<\n\n";
#endif
		/** @todo parsing molecule information */
		///////////////////////////////////////////////////////////////////////////
		MoleculeType& tp = molecules.new_type();
		tp.set_name(parse_attr<simString>(n, "name"));
		tp.set_idx_def(parse_attr<simI1>(n, "id"));
		tp.set_eva_num(parse_attr<simI1>(n, "amount"));

		///////////////////////////////////////////////////////////////////////////
		// molecule: geometry
		// using lambda function for parse molecule shape
		simVI3 rpos; // relative dots positions
		auto rpos_func = [&](node_t* sn /* node 'rpos' */) {
			simI1 i = parse_attr<simI1>(sn, "id");
			simI2 v = parse_str <simI2>(sn->value());
			simI3 r;
			r.x = v.x;
			r.y = v.y;
			r.z = i;
			rpos.push_back(r);
		};
		for_each_node(n->first_node("shape"), "rpos", rpos_func);
		// assign and clean up
		tp.set_dot_pos(rpos);
		//rpos.clear();

		///////////////////////////////////////////////////////////////////////////
		// molecule: bonding
		//MoleculeType::simBonds bonds;
		//auto bond_func = [&](node_t* sn /* node 'bond' */) {
		//	bonds.emplace_back();
		//	MoleculeType::simOneBond& new_bond = bonds.back();
		//	new_bond.energy = parse_attr<simF1>(sn, "energy");
		//	new_bond.target = parse_attr<simI1>(sn, "target");
		//	// get all bond rpos
		//	auto bond_pos_func = [&](node_t* nrpos /* node 'rpos' */) {
		//		simI2 v = parse_str <simI2>(nrpos->value());
		//		simI1 i = parse_attr<simI1>(nrpos, "id");
		//		new_bond.rpos.emplace_back();
		//		new_bond.rpos.back().x = v.x;
		//		new_bond.rpos.back().y = v.y;
		//		new_bond.rpos.back().z = i;
		//	};
		//	for_each_node(sn, "rpos", bond_pos_func);
		//};
		//for_each_node(n, "bond", bond_func);
		//// assign and clean up
		//tp.set_bond(bonds);
		//bonds.clear();
		///////////////////////////////////////////////////////////////////////////
	}

	/////////////////////////////////////////////////////////////////////////////
	// initialize substrate
	void init_substrate(const reader::node_t* n)
	{
		simI2 v = reader::parse_str<simI2>(n->value());
#ifndef NDEBUG
		cout << " initialize substrate: ";
		cout << "(" << v.x << "," << v.y << ")\n\n";
#endif
		sub.init(v.x, v.y);
	}
};

///////////////////////////////////////////////////////////////////////////////
// initialization function
void simula::init(const simChar* input)
{
	using namespace reader;
	// we cannot initialize a document with ptr, so open file here
	doc_t doc;
	xml_t xml(input);
	doc.parse<0>(xml.data());
	make_doc(doc);

	// initialize user defined constants
	for_each_node(doc.first_node(), "molecule", init_molecule_type);
	for_each_node(doc.first_node(), "substrate", init_substrate, 1);
#ifndef NDEBUG
	for (simI1 i = 0; i < molecules.type_num(); ++i)
	{
		molecules.type(i).debug();
	}
#endif
	// initialize the simulation system
	// initialize molecules
	for (simI1 i = 0; i < molecules.type_num(); ++i)
		init_molecule(molecules.type(i));
	// initialize random seed
	initRandSeed();
}
