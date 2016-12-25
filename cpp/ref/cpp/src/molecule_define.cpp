#include "molecule/molecule_define.h"

using namespace simula;

// here we land molecule index * 100 + dot indices on substrate */
const simSize Molecule::land(simSize i) const
{
	simSize id = 0;
	if (i < m_type->dot_num())
	{
		id = self_id() * max_dot_id + m_type->dot_pos()[i].z;
	}
#ifndef NDEBUG
	else
	{
		cerr << "ERROR: landing index out of bound" << endl;
	}
#endif
	return id;
}

#ifndef NDEBUG
/** @brief debug function */
void Molecule::debug()
{
	cout << " * molecule:";
	cout << " self_id: " << m_core.self_id;
	cout << " type_id: " << m_core.type_id;
	cout << " (" << m_core.x << "," << m_core.y << ")\n";
}
#endif
