#include <cstdlib>

#include "init.h"
#include "molecule.h"
#include "substrate.h"

#include "kernel.h"

using namespace simula;

int main(int argc, char *argv[])
{
	if (argc != 2) {
#ifndef NDEBUG
		// insuffcient argument number
		std::cerr << "the program need exactly one argument" << std::endl;
#endif
		return EXIT_FAILURE;
	}

	// normal argument number
	init(argv[1]);
#ifndef NDEBUG
	cout << "successfully initialized program" << endl;
#endif
#ifndef NDEBUG
	for (simSize id = 0; id < molecules.molecule_num(); ++id) {
		molecules.molecule(id).debug();
	}
#endif
#ifndef NDEBUG
	cout << "number of molecule types: " << molecules.type_num() << endl;
	cout << "number of molecule types: " << molecules.molecule_num() << endl;
#endif
	// land all points
	simI1 id = 0;
	while (id < molecules.molecule_num()) {
		simI1 Xpos = randInt(1, subXsize);
		simI1 Ypos = randInt(1, subYsize);
		if (sub.land(molecules.molecule(id), Xpos, Ypos, 0)) {
#ifndef NDEBUG
			molecules.molecule(id).debug();
#endif
			id++;
		}
	}
#ifndef NDEBUG
	cout << "successfully landed all molecules on thr substrate" << endl;
#endif
	// output substrate to file
	sub.print("output.txt");
	cout << sub << endl;
#ifndef NDEBUG
	cout << "successfully printed substrate" << endl;
#endif
	// @todo hopping


	simCuda::main_temp();


	// @todo hopping
	return 0;
}
