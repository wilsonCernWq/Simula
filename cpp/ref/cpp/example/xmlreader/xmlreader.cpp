#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "../../utils/rapidxml_1_13/rapidxml/rapidxml.hpp"

using namespace rapidxml;
using namespace std;

void allattr(xml_node<>* Node) {
  for (xml_attribute<>* attr = Node->first_attribute();
       attr; attr = attr->next_attribute()) 
    {
      cout << "Attribute: " << attr->name() << "\n";
      cout << "---> Value: " << attr->value() << "\n";
    }
}

int main (int argc, char *argv[]) 
{
  
  if (argc != 2) {

    // error handling
    std::cerr << "the program need exactly one argument" << std::endl;
    return -1;

  } else {

    // read file into string
    std::ifstream infile(argv[1]);
    std::stringstream buffer;
    buffer << infile.rdbuf();
    infile.close();
    std::string text(buffer.str());

    // read xml using rapidxml
    xml_document<> doc;
    doc.parse<0>(&text[0]);

    // print out data
    xml_node<>* rNode = doc.first_node();
    
    // for all sub nodes
    for (xml_node<>* iNode = rNode->first_node(); 
	 iNode; iNode = iNode->next_sibling()) 
      {

	// show all attribute
	allattr(iNode);
	cout <<  "Node value: " << iNode->value() << " ";
	cout <<  "Node Name: " << iNode->name() << "\n";
	
	// Interate over the sub nodes
	for(xml_node<> * jNode = iNode->first_node(); 
	    jNode; jNode = jNode->next_sibling())
	  {

	    allattr(jNode);
	    cout <<  "Node value: " << jNode->value() << " ";
	    cout <<  "Node name: " << jNode->name() << "\n";
	  }
      
      }    

    return 0;

  }

}

