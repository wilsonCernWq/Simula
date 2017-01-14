#pragma once

#include <iostream>

#include "simulaBaseType.h"

namespace simula {

	/**
	 * This is the super class for all xml nodes
	 */
	class Node {
	public:
		virtual ~Node() {}
		virtual void attribute(CHAR str, CHAR content) {
			std::cerr << "Warning: undefined attribute <<" << str << ">> found\n";
		}
		virtual Node* value(CHAR str) {
			std::cerr << "Warning: undefined node <<" << str << ">> found\n";
			return nullptr;
		}
	};

	/**
	* This is the implementation class tempolate for all xml nodes
	*/
	template <typename Core>
	class NodeImplementation : public Node {
	protected:
		static std::vector<Core> list;
		Core* core = nullptr;
	public: 
		NodeImplementation() {
			list.emplace_back();
			core = &list.back();
		}
		~NodeImplementation() {
			if (core == nullptr) { delete core; }
		}
		virtual void attribute(CHAR str, CHAR content) { 
			Node::attribute(str, content); 
		}
		virtual Node* value(CHAR str) {
			return Node::value(str);
		}
	};
	template <typename Core>
	std::vector<Core> NodeImplementation<Core>::list = std::vector<Core>();

}