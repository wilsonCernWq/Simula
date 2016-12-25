///////////////////////////////////////////////////////////////////////////////
//
// CUDA kernel function definition
//
///////////////////////////////////////////////////////////////////////////////
#ifndef _SIMULA_KERNEL_
#define _SIMULA_KERNEL_

#include "substrate.h"

///////////////////////////////////////////////////////////////////////////////
// project namespace
namespace simula {
	/////////////////////////////////////////////////////////////////////////////
	// namespace for CUDA kernel
	namespace simCuda {
		///////////////////////////////////////////////////////////////////////////
		// type for configuration (bond/shape)
		struct kRpos { simI1 x, y, i; };
		///////////////////////////////////////////////////////////////////////////
		// type for defining one bond
		struct kBond {
			simF1 energy; // energy
			simI1 target; // molecule target
			simI1   type; // bond identity index
			simSize size; // number of rpos for the bond 
		};
		///////////////////////////////////////////////////////////////////////////
		//
		// corresponding to Molecule
		//
		///////////////////////////////////////////////////////////////////////////
		struct kMolecule {
			simI1 x, y, d, i, t;
		};
		
		int main_temp();

	};
};


#endif // _SIMULA_KERNEL_