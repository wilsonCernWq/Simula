#pragma once

#include <iostream>
#include <memory>

#include "simulaBaseType.h"

namespace simula {
	/**
	 * This is the super class for all xml nodes
	 */
	class Node {
	public:
		virtual ~Node() {}

		virtual void attribute(CHAR str, CHAR content)
		{
			std::cerr << "Warning: undefined attribute <<" << str << ">> found\n";
		}

		virtual std::shared_ptr<Node> value(CHAR str)
		{
			std::cerr << "Warning: undefined node <<" << str << ">> found\n";
			return nullptr;
		}

	};
}