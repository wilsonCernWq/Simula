#include "molecule/molecule_type_define.h"

using namespace simula;

/** @brief set relative component position with a list of integer pairs **/
void MoleculeType::set_dot_pos(const simVI3& v)
{
	// set size and allocate array
	m_core.dot_num = v.size();
	m_core.dot_pos = new simI3[v.size()];
	// deep copying
	for (simSize i = 0; i < v.size(); ++i) { m_core.dot_pos[i] = v[i]; }
	return;
}
void MoleculeType::set_dot_pos(const simSize size, const simI3* value)
{
	// set size and allocate array
	m_core.dot_num = size;
	m_core.dot_pos = new simI3[size];
	// deep copying
	for (simSize i = 0; i < size; ++i) {	m_core.dot_pos[i] = value[i];	}
}

#ifndef NDEBUG
/** @brief  properties **/
void MoleculeType::debug()
{
	cout << "\n==> Molecule\n";
	cout << "   @ name: "    << m_name         << endl;
	cout << "   @ eva_num: " << m_core.eva_num << endl;
	cout << "   @ idx_def: " << m_core.idx_def << endl;
	cout << "   @ idx_gen: " << m_core.idx_gen << endl;
	cout << "   @ idx_num: " << m_core.dot_num << endl;
	cout << "   @ dot pos:\n";
	for (simSize i = 0; i < this->dot_num(); ++i) 
	{
		cout << "    # " << m_core.dot_pos[i].z;
		cout << " (" << m_core.dot_pos[i].x << "," << m_core.dot_pos[i].y << ")\n";
	}
	cout << endl;
}
#endif
