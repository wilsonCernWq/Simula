#include "global/global_function.h"

#include <cstdlib> // srand, rand
#include <ctime>   // time

using namespace simula;

/** @brief initialize random seed */
void simula::initRandSeed()
{
	srand(time(NULL));
}

/**
 * @brief generate integer random number within [lowerBound, upperBound)
 * @param lowerBound
 * @param upperBound
 */
simI1 simula::randInt(simI1 lowerBound, simI1 upperBound)
{
	return static_cast<simI1>(rand()) % (upperBound - lowerBound) + lowerBound;
}
