#pragma once

#include <engine/util/public/std.hpp>
#include <glm/vec4.hpp>
#include <glm/vec3.hpp>
#include <glm/vec2.hpp>
#include <glm/mat4x4.hpp>

RT_BEGIN_NAMESPACE

using vec2_t = glm::vec2;
using vec3_t = glm::vec3;
using vec4_t = glm::vec4;
using mat3_t = glm::mat3;
using mat4_t = glm::mat4;

enum collision_shape_type_e : u8_t
{
	COLLISION_SHAPE_TYPE_SPHERE,
	COLLISION_SHAPE_TYPE_PLANE,
};

struct ray_t
{
	vec3_t origin;
	vec3_t direction;
};

struct ray_hit_t
{
	vec3_t point;
	vec3_t normal;
	f32_t t;
};

struct sphere_t
{
	vec3_t center;
	f32_t radius;
};

struct plane_t
{
	vec3_t normal;
	f32_t offset;
};

RT_END_NAMESPACE