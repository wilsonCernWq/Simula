///////////////////////////////////////////////////////////////////////////////
//
// molecule list & type list definition and their handling functions
//
///////////////////////////////////////////////////////////////////////////////
#ifndef _SIMULA_MOLECULE_
#define _SIMULA_MOLECULE_

#include "molecule/molecule_define.h"

///////////////////////////////////////////////////////////////////////////////
// project namespace
namespace simula {
	/////////////////////////////////////////////////////////////////////////////
	// List CLass
	class MoleculeList {
	private:
		::std::vector<Molecule>     m_main_list;
		::std::vector<MoleculeType> m_type_list;
	public:

		~MoleculeList() { clean(); }

		// cleaning heap memory
		inline void clean() {
			for (simSize i = 0; i < m_type_list.size(); ++i) { m_type_list[i].clean(); }
			for (simSize i = 0; i < m_main_list.size(); ++i) { m_main_list[i].clean(); }
		}

		// generate a new MoleculeType/Molecule
		MoleculeType& new_type();
		Molecule& new_molecule(const MoleculeType& type);

		// get the MoleculeType by index
		inline MoleculeType& type(simSize i) { return m_type_list[i]; }

		// get the Molecule by index
		inline Molecule& molecule(simSize i) { return m_main_list[i]; }

		// get the number of MoleculeType/Molecule
		inline const simSize type_num() { return m_type_list.size(); }
		inline const simSize molecule_num()	{ return m_main_list.size(); }

	};

	/////////////////////////////////////////////////////////////////////////////
	// global variable
	extern MoleculeList molecules;
};

#endif // _SIMULA_MOLECULE_
