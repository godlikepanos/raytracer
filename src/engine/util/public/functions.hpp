#pragma once

#include "std.hpp"

RT_BEGIN_NAMESPACE

inline f32_t rand_0f_1f()
{
	return (f32_t)drand48();
}

template<typename t>
t* new_array(u32_t size)
{
	assert(size > 0);
	return (t*)malloc(sizeof(t) * size);
}

RT_END_NAMESPACE