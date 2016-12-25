#ifndef _SIMULA_MOLECULE_DEFINE_
#define _SIMULA_MOLECULE_DEFINE_

#include "molecule_type_define.h"

///////////////////////////////////////////////////////////////////////////////
// project namespace
namespace simula {

	/////////////////////////////////////////////////////////////////////////////
	// Molecule class defining every simulation object
	//  var: type
	//  var: type_id
	//	var: self_id
	//	var: x
	//	var: y
	//	var: d
	class Molecule {
	public:
		static const simSize max_dot_id = 100;
		/////////////////////////////////////////////////////////////////////////////
		// convert land id to molecule self id
		static inline simI1 land_to_sid(simI1 v)
		{
			return static_cast<simI1>(v / max_dot_id);
		}
		/////////////////////////////////////////////////////////////////////////////
		// @brief convert land id to molecule rpos id
		static inline simI1 land_to_rid(simI1 v)
		{
			return static_cast<simI1>(v % max_dot_id);
		}
	public:
		struct core_t {
			simSize self_id; // index of the instance
			simSize type_id; // index of the type
			simI1 d; // molecular direction
			simI1 x; // molecular x-position
			simI1 y; // molecular y-position
		};
	private:
		const MoleculeType * m_type = NULL; // pointer to its molecule type
		Molecule::core_t     m_core;
	public:

		// cleaning heap memory
		inline void clean() {}

		///////////////////////////////////////////////////////////////////////////////
		// Getter
		inline const simSize type_id() const { return m_core.type_id; }
		inline const simSize self_id() const { return m_core.self_id; }
		inline const simI1 x() const { return m_core.x; }
		inline const simI1 y() const { return m_core.y; }
		inline const simI1 d() const { return m_core.d; }
		inline const MoleculeType& type() const { return *m_type; }
		const simSize land(simSize i) const;

		///////////////////////////////////////////////////////////////////////////////
		// Setter
		inline void set_type(const MoleculeType& p) { m_type = &p; }
		inline void set_type_id(simSize v) { m_core.type_id = v; }
		inline void set_self_id(simSize v) { m_core.self_id = v; }
		inline void set_x(const simI1 x) { m_core.x = x; }
		inline void set_y(const simI1 y) { m_core.y = y; }
		inline void set_d(const simI1 d) { m_core.d = d; }

#ifndef NDEBUG
		// debug function
		void debug();
#endif

	};
};

#endif // _SIMULA_MOLECULE_DEFINE_
