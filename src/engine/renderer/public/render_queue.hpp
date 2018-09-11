#pragma once

#include "common.hpp"
#include "material.hpp"
#include <engine/math/public/types.hpp>

RT_BEGIN_NAMESPACE

struct renderable_t
{
	union
	{
		collision_shape_t base;
		sphere_t sphere;
		plane_t plane;
	};

	material_t material;
};

struct render_queue_t
{
	renderable_t* renderables;
	u32_t renderable_count;
};

RT_END_NAMESPACE
