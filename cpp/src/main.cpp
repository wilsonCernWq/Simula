#include <vector>

typedef int          mInt ;
typedef double       mReal;
typedef unsigned int mSize;

struct mTarget {
	mInt rx, ry, rz, rd, si, sj;
};

struct mReaction {
	int gIdx, //< global index 
		rIdx; //< relative (local) index
	std::vector<mTarget> targets;
};

struct mComponent {
	int idx;
	int rx, ry, rz, rd;
	std::vector<int> directions;
	std::vector<mReaction> reactions;
};

class molecule {
public:
	int x, y, z, d;
	std::vector<mComponent> components;
};

// Matrices are stored in row-major order:
// M(row, col) = *(M.elements + row * M.width + col)
struct substrate {
	int width;
	int height;
	float* elements;
};


int main() {

	return 0;
};