#pragma once

#include "common.hpp"

RT_BEGIN_NAMESPACE

enum collision_shape_type : u8
{
	SPHERE,
	PLANE,
};

struct collision_shape
{
	collision_shape_type type;
};

struct ray
{
	vec3 origin;
	vec3 direction;
};

struct ray_hit
{
	vec3 point;
	vec3 normal;
	f32 t;
};

struct sphere
{
	collision_shape base;
	vec3 center;
	f32 radius;
};

struct plane
{
	collision_shape base;
	vec3 normal;
	f32 offset;
};

RT_END_NAMESPACE