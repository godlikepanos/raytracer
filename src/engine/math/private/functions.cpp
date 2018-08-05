#include <engine/math/public/functions.hpp>
#include <engine/math/public/types.hpp>

RT_BEGIN_NAMESPACE

bool ray_cast_sphere(const ray& ray, const sphere& sphere, f32 t_min, f32 t_max, ray_hit* hit)
{
	vec3 oc = ray.origin - sphere.center;
	f32 a = dot(ray.direction, ray.direction);
	f32 b = dot(oc, ray.direction);
	f32 c = dot(oc, oc) - sphere.radius * sphere.radius;
	f32 d = b * b - a * c;

	if(d > 0.0f)
	{
		f32 tmp = (-b - sqrt(d)) / a;
		if(tmp < t_max && tmp > t_min)
		{
			hit->t = tmp;
			hit->point = ray.origin + ray.direction * hit->t;
			hit->normal = (hit->point - sphere.center) / sphere.radius;
			return true;
		}
	}

	return false;
}

bool ray_cast_plane(const ray& ray, const plane& plane, f32 t_min, f32 t_max, ray_hit* hit)
{
	f32 d = plane_point_distance(plane, ray.origin);
	f32 a = dot(plane.normal, ray.direction);

	if(d > 0.0f && a < 0.0f)
	{
		f32 tmp = -d / a;
		if(tmp < t_max && tmp > t_min)
		{
			hit->t = -d / a;
			assert(hit->t > 0.0f);
			hit->point = ray.origin + ray.direction * hit->t;
			hit->normal = plane.normal;
			return true;
		}
	}

	return false;
}

bool ray_cast(const ray& ray, const collision_shape& shape, f32 t_min, f32 t_max, ray_hit* hit)
{
	bool intersects;
	switch(shape.type)
	{
	case collision_shape_type::SPHERE:
		intersects = ray_cast_sphere(ray, *(const sphere*)&shape, t_min, t_max, hit);
		break;
	case collision_shape_type::PLANE:
		intersects = ray_cast_plane(ray, *(const plane*)&shape, t_min, t_max, hit);
		break;
	default:
		assert(0);
		intersects = false;
	}

	return intersects;
}

vec3 rand_direction_in_cone(const vec3& cone_dir, f32 cone_angle)
{
	f32 phi = (f32)rand() / RAND_MAX; // angle to z
	phi = phi * 2.0f - 1.0f;
	phi *= cone_angle;
	phi = M_PI / 2.0f - phi;

	f32 theta = (f32)rand() / RAND_MAX; // angle to x
	theta = theta * 2.0f - 1.0f;
	theta *= cone_angle;

	vec3 rand_dir; // around x axis
	rand_dir.x = cos(theta) * sin(phi);
	rand_dir.y = sin(theta) * sin(phi);
	rand_dir.z = cos(phi);

	// Create a rot matrix
	vec3 x_axis = cone_dir;
	vec3 y_axis = vec3(0.0f, 1.0f, 0.0f);
	vec3 z_axis = cross(x_axis, y_axis);
	y_axis = cross(z_axis, x_axis);

	mat3 rot(x_axis, y_axis, z_axis);

	return rot * rand_dir;
}

RT_END_NAMESPACE