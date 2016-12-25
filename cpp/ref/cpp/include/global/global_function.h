///////////////////////////////////////////////////////////////////////////////
//
// global function definition
//
///////////////////////////////////////////////////////////////////////////////
#ifndef _SIMULA_GLOBAL_FUNCTION_
#define _SIMULA_GLOBAL_FUNCTION_

#include "global/global_type.h"

namespace simula {

	/** @brief initialize random seed */
	void initRandSeed();

	/**
	 * @brief generate integer random number within [lowerBound, upperBound)
	 * @param lowerBound
	 * @param upperBound
	 */
	simI1 randInt(simI1 lowerBound, simI1 upperBound);

	/**
	 * @brief positive modulo
	 * @param x Input value
	 * @param n Base integer
	 **/
	inline simI1 pmod(simI1 x, simI1 n)
	{
		return (x % n + n) % n;
	}

};

#endif // _SIMULA_GLOBAL_FUNCTION_
