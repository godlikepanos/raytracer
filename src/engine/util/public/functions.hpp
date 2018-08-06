#pragma once

#include "std.hpp"

RT_BEGIN_NAMESPACE

f32 rand_0f_1f()
{
	return (f32)rand() / RAND_MAX;
}

RT_END_NAMESPACE