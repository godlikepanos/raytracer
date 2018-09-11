#pragma once

#include "types.hpp"
#include <glm/geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>

RT_BEGIN_NAMESPACE

using glm::dot;
using glm::clamp;
using glm::perspective;

inline void sphere_init(sphere_t& s)
{
	s.base.type = collision_shape_type_e::SPHERE;
}

inline void plane_init(plane_t& p)
{
	p.base.type = collision_shape_type_e::PLANE;
}

inline f32_t plane_point_distance(const plane_t& plane, const vec3_t& point)
{
	return dot(plane.normal, point) - plane.offset;
}

// Ray intersection
bool_t ray_cast_sphere(const ray_t& ray, const sphere_t& sphere, f32_t t_min, f32_t t_max, ray_hit_t& hit);
bool_t ray_cast_plane(const ray_t& ray, const plane_t& plane, f32_t t_min, f32_t t_max, ray_hit_t& hit);

bool_t ray_cast(const ray_t& ray, const collision_shape_t& shape, f32_t t_min, f32_t t_max, ray_hit_t& hit);

// Other
vec3_t rand_direction_in_cone(const vec3_t& cone_dir, f32_t cone_angle);
vec3_t rand_point_in_unit_sphere();

RT_END_NAMESPACE
