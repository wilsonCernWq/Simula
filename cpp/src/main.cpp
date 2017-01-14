#include <iostream>
#include <fstream>
#include <sstream>

#include "xmlReader.h"
#include "xmlNodeDefinition_Action.h"

/** Main program to parse the input XML file
 */
int main(int argc, char *argv[]) {

	if (argc != 2) { ///< error handling (Too many/few arguments)

		std::cerr << "the program need exactly one argument" << std::endl;
		return -1;

	} else {

		/// parse file into string
		std::ifstream infile(argv[1]);
		std::stringstream buffer;
		buffer << infile.rdbuf();
		infile.close();

		/// parse xml file
		simula::parseInput(buffer.str());

		return 0;

	}

}
