#ifndef _SIMULA_SUBSTRATE_
#define _SIMULA_SUBSTRATE_

#include "molecule.h"

///////////////////////////////////////////////////////////////////////////////
// local variable namespace
namespace m_substrate {
	extern const simula::simI1 m_bg; //< background value
};

///////////////////////////////////////////////////////////////////////////////
// project namespace
namespace simula {
	/////////////////////////////////////////////////////////////////////////////
	// Substrate class defining substrate behaviors
	class Substrate {
	private:
		simBool m_init_flag = false;
		simSize m_xlen; //< X dimension size
		simSize m_ylen; //< Y dimemsion size
		simVI1  m_data; //< array of pointers to Molecule
	private:

	public:

		///////////////////////////////////////////////////////////////////////////
		// Getter
		inline simI1* data() { return m_data.data(); }
		inline simI1  xlen() { return m_xlen; }
		inline simI1  ylen() { return m_ylen; }

		///////////////////////////////////////////////////////////////////////////
		// Initialize substrate with value background value
		void init(simSize Xsize, simSize Ysize);

		///////////////////////////////////////////////////////////////////////////
		// Point value getter
		const simI1 get_sub(const simI1 x, const simI1 y) const;

		///////////////////////////////////////////////////////////////////////////
		// Point value setter
		void set_sub(const simI1 x, const simI1 y, const simI1 value);

		///////////////////////////////////////////////////////////////////////////
		// Check if the point is empty
		inline simBool is_empty(const simI1 x, const simI1 y) const
		{
			return (get_sub(x, y) == m_substrate::m_bg);
		}

		///////////////////////////////////////////////////////////////////////////
		// Check if the relative positions are all empty
		simBool is_empty(const MoleculeType& m, const simI1 xc, const simI1 yc) const;

		///////////////////////////////////////////////////////////////////////////
		// Land molecule on the position
		simBool land(Molecule& m, const simI1 xc, const simI1 yc, const simI1 dc);

		///////////////////////////////////////////////////////////////////////////
		// print substrate into file
		void print(const simChar* name);

		///////////////////////////////////////////////////////////////////////////
		// overload output function
		friend std::ostream& operator<<(std::ostream& os, const Substrate& sub);

	};

	///////////////////////////////////////////////////////////////////////////
	// overload output function
	std::ostream& operator<<(std::ostream& os, const Substrate& sub);

	///////////////////////////////////////////////////////////////////////////
	// local variable
	extern simula::Substrate sub;
};

#endif // _SIMULA_SUBSTRATE_
