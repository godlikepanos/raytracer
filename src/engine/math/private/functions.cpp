#include <engine/math/public/functions.hpp>
#include <engine/math/public/types.hpp>
#include <engine/util/public/functions.hpp>
#include <cstdio>

RT_BEGIN_NAMESPACE

bool_t ray_cast_sphere(const ray_t& ray, const sphere_t& sphere, f32_t t_min, f32_t t_max, ray_hit_t& hit)
{
	vec3_t oc = ray.origin - sphere.center;
	f32_t a = dot(ray.direction, ray.direction);
	f32_t b = dot(oc, ray.direction);
	f32_t c = dot(oc, oc) - sphere.radius * sphere.radius;
	f32_t d = b * b - a * c;

	if(d > 0.0f)
	{
		f32_t tmp = (-b - sqrt(d)) / a;
		if(tmp < t_max && tmp > t_min)
		{
			hit.t = tmp;
			hit.point = ray.origin + ray.direction * hit.t;
			hit.normal = (hit.point - sphere.center) / sphere.radius;
			return true;
		}
	}

	return false;
}

bool_t ray_cast_plane(const ray_t& ray, const plane_t& plane, f32_t t_min, f32_t t_max, ray_hit_t& hit)
{
	f32_t d = plane_point_distance(plane, ray.origin);
	f32_t a = dot(plane.normal, ray.direction);

	if(d > 0.0f && a < 0.0f)
	{
		f32_t tmp = -d / a;
		if(tmp < t_max && tmp > t_min)
		{
			hit.t = -d / a;
			assert(hit.t > 0.0f);
			hit.point = ray.origin + ray.direction * hit.t;
			hit.normal = plane.normal;
			return true;
		}
	}

	return false;
}

bool_t ray_cast(const ray_t& ray, const collision_shape_t& shape, f32_t t_min, f32_t t_max, ray_hit_t& hit)
{
	bool_t intersects;
	switch(shape.type)
	{
	case collision_shape_type_e::SPHERE:
		intersects = ray_cast_sphere(ray, *(const sphere_t*)&shape, t_min, t_max, hit);
		break;
	case collision_shape_type_e::PLANE:
		intersects = ray_cast_plane(ray, *(const plane_t*)&shape, t_min, t_max, hit);
		break;
	default:
		assert(0);
		intersects = false;
	}

	return intersects;
}

vec3_t rand_direction_in_cone(const vec3_t& cone_dir, f32_t cone_angle)
{
	// https://math.stackexchange.com/questions/56784/generate-a-random-direction-within-a-cone/205589#205589

	// Create a random vector where the cone angle is in +z
	f32_t cos_cone_angle = cos(cone_angle);
	f32_t z = rand_0f_1f() * (1.0f - cos_cone_angle) + cos_cone_angle;
	f32_t phi = rand_0f_1f() * 2.0f * M_PI;
	f32_t sq = sqrt(1.0f - z * z);
	f32_t x = sq * cos(phi);
	f32_t y = sq * sin(phi);

	// Create a rot matrix
	vec3_t z_axis = cone_dir;
	vec3_t x_axis = vec3_t(1.0f, 0.0f, 0.0f);
	vec3_t y_axis = cross(z_axis, x_axis);
	x_axis = cross(y_axis, z_axis);

	mat3_t rot(x_axis, y_axis, z_axis);

	return rot * vec3_t(x, y, z);
}

vec3_t rand_point_in_unit_sphere()
{
	vec3_t p;
	do
		p = 2.0f * vec3_t(rand_0f_1f(), rand_0f_1f(), rand_0f_1f()) - vec3_t(1.0f);
	while(p.length() <= 1.0f);

	return p;
}

RT_END_NAMESPACE
