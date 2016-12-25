///////////////////////////////////////////////////////////////////////////////
//
// substrate class definition
//
///////////////////////////////////////////////////////////////////////////////
#include "substrate.h"

#include <ostream> // ostream
#include <fstream> // ofstream

///////////////////////////////////////////////////////////////////////////////
// used namespace
using namespace m_substrate;
using namespace simula;

///////////////////////////////////////////////////////////////////////////////
// initialize variables from header file
const simI1 m_substrate::m_bg = 0; // background value
Substrate simula::sub;

///////////////////////////////////////////////////////////////////////////////
// Initialize substrate with value background value
void Substrate::init(simSize Xsize, simSize Ysize)
{
	if (!m_init_flag)
	{
		m_init_flag = true;
		m_xlen = Xsize;
		m_ylen = Ysize;
		for (simI1 i = 0; i < m_xlen * m_ylen; ++i)
		{
			m_data.push_back(m_bg);
		}
#ifndef NDEBUG
		cout << " **** generating substrate (";
		cout << m_xlen << "," << m_ylen << ") ****\n\n";
#endif
		// need to be modified once the substrate class is fully implemented
		subXsize = m_xlen;
		subYsize = m_ylen;
	}
#ifndef NDEBUG
	else
	{
		cerr << " ERROR: substrate is already generated\n\n";
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////
// Point value getter
const simI1 Substrate::get_sub(const simI1 x, const simI1 y) const
{
	simI1 mx = pmod(x, m_xlen);
	simI1 my = pmod(y, m_ylen);
	return m_data[my * m_xlen + mx];
}

///////////////////////////////////////////////////////////////////////////////
// Point value setter
void Substrate::set_sub(const simI1 x, const simI1 y, const simI1 value)
{
	simI1 mx = pmod(x, m_xlen);
	simI1 my = pmod(y, m_ylen);
	m_data[my * m_xlen + mx] = value;
}

///////////////////////////////////////////////////////////////////////////////
// Check if the relative positions are all empty
simBool Substrate::is_empty
(const MoleculeType& tp, const simI1 xc, const simI1 yc) const
{
	simBool empty = true;
	for (simI1 i = 0; i < tp.dot_num(); ++i) {
		simI1 x = xc + tp.dot_pos()[i].x;
		simI1 y = yc + tp.dot_pos()[i].y;
		if (!is_empty(x, y)) { empty = false; break; }
	}
	return empty;
}

///////////////////////////////////////////////////////////////////////////////
// Land molecule on the position
simBool Substrate::land(Molecule& m, const simI1 xc, const simI1 yc, const simI1 dc)
{
	// check overlap
	if (is_empty(m.type(), xc, yc))
	{
		// set molecular position & direction
		m.set_x(xc);
		m.set_y(yc);
		m.set_d(dc);
		// set point values on substrate
		for (simI1 i = 0; i < m.type().dot_num(); ++i) {
			simI1 x = xc + m.type().dot_pos()[i].x;
			simI1 y = yc + m.type().dot_pos()[i].y;
			set_sub(x, y, m.land(i) /* land index containing molecule index and rpos index */);
		}
		return true;
	}
	else {
#ifndef NDEBUG
		cout << "landing fail" << endl;
#endif
		// make no changes for checking failure
		return false;
	}
}

///////////////////////////////////////////////////////////////////////////////
// print substrate into file
void Substrate::print(const simChar* name)
{
	std::ofstream file(name);
	if (file.is_open()) 
	{
		file << *this;
		file.close();
	}
	else {
#ifndef NDEBUG
		std::cerr << "fail to open a file" << std::endl;
#endif
	}
}

///////////////////////////////////////////////////////////////////////////////
// overload output function
std::ostream& simula::operator<<(std::ostream& os, const Substrate& sub)
{
	for (simI1 i = 0; i < sub.m_xlen; ++i) {
		for (simI1 j = 0; j < sub.m_ylen; ++j) {
			// check if the point is occupied by something
			simI1 v = Molecule::land_to_sid(sub.get_sub(i,j));
			// if there is nothing, filled by zero
			if (v == m_bg) {
				os << 0 << " ";
			}
			// if there is something, fill its type index instead the self_id (reference index/data index)
			else {
				// since all indices are counted from 1, we need to substract 1 first
				os << molecules.molecule(v - 1).type_id() << " ";
			}
		}
		os << std::endl;
	}
	return os;
}

