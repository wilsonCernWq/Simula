#pragma once

#include <iostream>

#include "simulaBaseType.h"

namespace simula {

	/** 
	 * \brief This is the super class for all xml nodes
	 */
	class Node {
	public:
		virtual ~Node() {}
		virtual void attribute(CHAR str, CHAR content) {
			std::cerr << "Warning: undefined attribute << " << str << " >> found\n";
		}
		virtual Node* value(CHAR str) {
			std::cerr << "Warning: undefined node << " << str << " >> found\n";
			return nullptr;
		}
	};

	/**
	* \brief This is the implementation class tempolate for all xml nodes
	* For every node implementation, one should first define the core data structure,
	* and then plug-in the core data structure into this template. 
	* The template class will define a static vector array for storing the core data,
	* this is particularly designed for GPU/CUDA acceleration. During the simulation 
	* process, the array containing only core structure will be extracted, and mapped 
	* into lower level C interface.
	*/
	template <typename Core>
	class NodeImplementation : public Node {
	protected:
		static std::vector<Core> list;
		Core* core = nullptr; ///< pointer to the core data
	public: 
		NodeImplementation() {
			list.emplace_back();
			core = &list.back();
		}
		~NodeImplementation() {
			if (core == nullptr) { delete core; }
		}
		///< Avoid namespace shadowing 
		virtual void attribute(CHAR str, CHAR content) { 
			Node::attribute(str, content); 
		}
		virtual Node* value(CHAR str) {
			return Node::value(str);
		}
	};
	///< Initialization for static vector field
	template <typename Core>
	std::vector<Core> NodeImplementation<Core>::list = std::vector<Core>();

}