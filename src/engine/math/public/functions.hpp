#pragma once

#include "ray.hpp"
#include "sphere.hpp"

RT_BEGIN_NAMESPACE

inline f32 ray_hit_sphere(const ray& ray, const sphere& sphere)
{
	vec3 oc = ray.origin - sphere.center;
	f32 a = dot(ray.direction, ray.direction);
	f32 b = 2.0f * dot(oc, ray.direction);
	f32 c = dot(oc, oc) - sphere.radius * sphere.radius;
	f32 d = b * b - 4.0f * a * c;
	if(d < 0.0f)
		return -1.0f;
	else
		return (-b - sqrt(d)) / (2.0f * a);
}

RT_END_NAMESPACE