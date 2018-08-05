#pragma once

#include "types.hpp"

RT_BEGIN_NAMESPACE

inline void sphere_init(sphere* s)
{
	s->base.type = collision_shape_type::SPHERE;
}

inline void plane_init(plane* p)
{
	p->base.type = collision_shape_type::PLANE;
}

inline f32 plane_point_distance(const plane& plane, const vec3& point)
{
	return dot(plane.normal, point) - plane.offset;
}

// Ray intersection
bool ray_cast_sphere(const ray& ray, const sphere& sphere, f32 t_min, f32 t_max, ray_hit* hit);
bool ray_cast_plane(const ray& ray, const plane& plane, f32 t_min, f32 t_max, ray_hit* hit);

bool ray_cast(const ray& ray, const collision_shape& shape, f32 t_min, f32 t_max, ray_hit* hit);

RT_END_NAMESPACE