#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <../../utils/rapidxml_1_13/rapidxml/rapidxml.hpp>

namespace simula {

  /// define internal types inside the program
  typedef std::string  CHAR;
  typedef int          INT;
  typedef unsigned int UINT;
  typedef double       REAL;

  /// define low level data structures (for CUDA)
  struct Action_Core {
    UINT gId; ///< global index (used in KMC lookup)
    UINT rId; ///< relative (local) index
    UINT id_component;
    INT state_i; ///< initial state
    INT state_f; ///< final state
  };
  std::vector<Action_Core> vector_actions;

  struct Occurrence_Core {
    UINT gIdx; ///< global index (used in KMC lookup)
    UINT rIdx; ///< relative (local) index
    UINT occur_min, occur_max;
    UINT iid_action, fid_action; ///< initial action index, final action index
  };
  std::vector<Occurrence_Core> vector_occurrences;
  
  /** Super class
   */
  class Node {
  public:
    virtual ~Node () {}

    virtual void attribute (CHAR str, CHAR content) 
    {
      std::cerr << "Warning: undefined attribute <<" << str << ">> found\n";
    }

    virtual Node* value (CHAR str) 
    { 
      std::cerr << "Warning: undefined node <<" << str << ">> found\n";
      return NULL;
    }

  };

  /**      
   */
  class Action : public Node {
  public:
    Action_Core* core = NULL;
  public:
    Action() {
      core = new Action_Core();
    }
    ~Action() {
      if (core != NULL) { delete core; }
    }

    void attribute (CHAR str, CHAR content) {
      if (str.compare("cid") == 0) { 
	core->id_component = stoi(content); 
      }
      else if (str.compare("initial") == 0) { 
        core->state_i = stoi(content); 
      }
      else if (str.compare("final") == 0)   { 
	core->state_f = stoi(content);
      }
    }

  };

  class Occurrence : public Node {
  public:
    Occurrence_Core* core = NULL;
    std::vector<Action*> actions;
  public:
    Occurrence() {
      core = new Occurrence_Core();
    }
    ~Occurrence() {
      if (core != NULL) { delete core; }
    }

    void attribute (CHAR str, CHAR content) {
      if (str.compare("min") == 0) { 
	core->min = stoi(content); 
      }
      else if (str.compare("max") == 0) { 
	core->max = stoi(content); 
      }
    }

    /// define new node
    /// 1) create new element inside the global list
    
    virtual Node* value (CHAR) {
      if (str.compare("action") == 0) { 
	auto action = new Action;
	this->actions.push_back(*action);
	return action;
      }
    }

  };

  class mSpot : public mNode {
  public:

    mIdx gIdx, ///< global index (used in KMC lookup)
      rIdx; ///< relative (local) index
    mStr cid;
    mInt rx, ry, rz, rd;

    /// filling attributes
    void attribute (std::string str, std::string content) {
      if (str.compare("cid") == 0) {
	this->cid = content;
      }
      else if (str.compare("rx") == 0) {
	this->rx = stoi(content);
      }
      else if (str.compare("ry") == 0) {
	this->ry = stoi(content);
      }
      else if (str.compare("rz") == 0) {
	this->rz = stoi(content);
      }
      else if (str.compare("rd") == 0) {
	this->rd = stoi(content);
      }
    }

  };

  class mReaction : public mNode {
  public:

    mIdx gIdx, ///< global index (used in KMC lookup)
      rIdx; ///< relative (local) index
    mReal energy; ///< reaction energ
    mInt initial, final;
    std::vector<mSpot> spots;
    std::vector<mOccurance> occurances;

    virtual mNode* value (std::string str) {
      if (str.compare("spot") == 0) { 
	auto node = new mSpot;
	this->spots.push_back(*node);
	return node;
      }
      else if (str.compare("occurance") == 0) { 
	auto node = new mOccurance;
	this->occurances.push_back(*node);
	return node;
      }
    }

  };

  /** Definition of a molecule's components.
   */
  class mComponent : public mNode {
  public:

    mIdx gIdx;
    mStr cIdx; ///< component index
    mInt bgcolor; ///< initial color (backgrounf color)
    std::vector<mReaction> reactions;

    virtual mNode* value (std::string str) {
      if (str.compare("reaction") == 0) { 
	auto node = new mReaction;
	this->reactions.push_back(*node);
	return node;
      }
    }

  };

  /** Definition of a molecule, the top-level structure   
   */
  class mMolecule : public mNode {
  public:

    int x, y, z, d;
    std::vector<mComponent> components;

  };

  /** Definition of a substrate
   * Matrices are stored in row-major order:
   * M(row, col) = *(M.elements + row * M.width + col)
   */
  class mSubstrate : public mNode {
  public:

    mInt width;
    mInt height;
    float* elements;

    /// filling attributes
    void attribute (std::string str, std::string content) {
      if (str.compare("width") == 0) {
	this->width = stoi(content);
      }
      else if (str.compare("height") == 0) {
	this->height = stoi(content);
      }
    }

  };

  class mDoc : public mNode {
  public:

    std::vector<mMolecule> molecules;
    std::vector<mComponent> components;
    mSubstrate* substrate = NULL;
  
    /// filling node values
    virtual mNode* value (std::string str) {

      if (str.compare("component") == 0) { 
	auto node = new mComponent;
	this->components.push_back(*node);
	return node;
      }
      else if (str.compare("molecule") == 0) { 
	auto node = new mMolecule;
	this->molecules.push_back(*node);
	return node;
      }
      else if (str.compare("substrate") == 0) {
	if (this->substrate == NULL) {
	  this->substrate = new mSubstrate;
	} else {
	  std::cerr << "Error: multiple substrate definitions found\n";
	}
	return substrate;
      }
    }

  };

  /** Function for reading XML input files
   */
  void readnode(rapidxml::xml_node<>* n, mInt level, mNode* parent) 
  {
    using namespace rapidxml;
    using namespace std;
    typedef xml_attribute<> xmlattr;
    typedef xml_node<> xmlnode;
    std::string indent = "   ";

    ///< node information
    for (mInt i = 0; i<level; ++i) { cout << indent; } ///< level indentation
    cout <<  "node: " << n->name() << " ";
    cout <<  "value: " << n->value() << "\n";

    ///< read attributes
    for (xmlattr* a = n->first_attribute(); a; a = a->next_attribute()) {
      for (mInt i = 0; i < level; ++i) { cout << indent; } ///< level indentation
      ///< format: [...] attribute: xxx value: xxx
      cout << "-attribute: " << a->name() << " ";
      cout << "value: " << a->value() << "\n";
    }

    /// TODO define objects here
    mNode* child = parent->value(string(n->name()));

    ///< read children
    for(xmlnode* c = n->first_node(); c; c = c->next_sibling()) {
      readnode(c,level+1,child);
    }

  }

}

/** Main program to parse the input XML file
 */
int main (int argc, char *argv[]) {
  
  if (argc != 2) { ///< error handling (Too many/few arguments)

    std::cerr << "the program need exactly one argument" << std::endl;
    return -1;

  } else {

    /// read file into string
    std::ifstream infile(argv[1]);
    std::stringstream buffer;
    buffer << infile.rdbuf();
    infile.close();
    std::string text(buffer.str());

    /// read xml using rapidxml
    rapidxml::xml_document<> doc;
    doc.parse<0>(&text[0]);

    /// print out data
    rapidxml::xml_node<>* rNode = doc.first_node();
    
    /// for all sub nodes
    mDoc* root = new mDoc;
    for(rapidxml::xml_node<>* c = rNode->first_node(); c; c = c->next_sibling()) {
      readnode(c,0,root);
    }

    return 0;

  }

}
