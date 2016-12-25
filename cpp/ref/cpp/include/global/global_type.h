///////////////////////////////////////////////////////////////////////////////
//
// define types used throughout the project
//
///////////////////////////////////////////////////////////////////////////////
#ifndef _SIMULA_GLOBAL_TYPE_
#define _SIMULA_GLOBAL_TYPE_

#include <boost/algorithm/string.hpp>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
// local variable namespace
namespace m_global {
	template<typename T1>
	struct T2
	{
		T1 x, y;
		T2() {}
		T2(T1 xi, T1 yi)
			: x(xi), y(yi) {}
	};
	template<typename T1>
	struct T3
	{
		T1 x, y, z;
		T3() {}
		T3(T1 xi, T1 yi, T1 zi)
			: x(xi), y(yi), z(zi) {}
	};
};

///////////////////////////////////////////////////////////////////////////////
// project namespace
namespace simula {

	/////////////////////////////////////////////////////////////////////////////
	// basic types
	typedef int    simI1;
	typedef double simF1;
	typedef bool   simBool;
	typedef char   simChar;
	/////////////////////////////////////////////////////////////////////////////
	// extended types
	typedef std::string simString;
	typedef size_t      simSize;

	/////////////////////////////////////////////////////////////////////////////
	// conbined types
	typedef ::m_global::T2<simI1> simI2;
	typedef ::m_global::T2<simF1> simF2;
	typedef ::m_global::T3<simI1> simI3;
	typedef ::m_global::T3<simF1> simF3;
	/////////////////////////////////////////////////////////////////////////////
	// list types
	typedef ::std::vector<simI1> simVI1;
	typedef ::std::vector<simI2> simVI2;
	typedef ::std::vector<simI3> simVI3;
	typedef ::std::vector<simF1> simVF1;
	typedef ::std::vector<simF2> simVF2;
	typedef ::std::vector<simF3> simVF3;
	typedef ::std::vector<simString> simStrVec;

};

#endif // _SIMULA_GLOBAL_TYPE_
