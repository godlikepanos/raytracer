#include <engine/misc/public/tga.hpp>
#include <engine/math/public/all.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <cstdio>

int main(int, char**)
{
	using namespace RT_NAMESPACE;

	const u32 width = 512;
	const u32 height = 256;

	vec3 cam_pos = vec3(0.0f, 1.0f, 0.0f);

	mat4x4 view_mat = lookAt(vec3(0.0f, 0.0f, -1.0f), cam_pos, vec3(0.0f, 1.0f, 0.0f));
	mat4x4 cam_trf = inverse(view_mat);
	mat4x4 proj_mat = perspective(glm::radians(45.0f), (f32)width / (f32)height, 0.1f, 100.0f);
	mat4x4 vp_mat = proj_mat * view_mat;
	mat4x4 inv_vp_mat = inverse(vp_mat);

	sphere sphere;
	sphere_init(&sphere);
	sphere.center = vec3(0.0f, 1.0f, -5.0f);
	sphere.radius = 1.0f;

	plane plane;
	plane_init(&plane);
	plane.offset = 0.0f;
	plane.normal = vec3(0.0f, 1.0f, 0.0f);

	u8 data[height][width][3];
	for(u32 w = 0; w < width; ++w)
		for(u32 h = 0; h < height; ++h)
		{
			vec2 ndc((w + 0.5f) / width, (h + 0.5f) / height);
			ndc = ndc * 2.0f - 1.0f;

			vec3 view_dir = inv_vp_mat * vec4(ndc, 1.0f, 1.0f);
			ray ray;
			ray.direction = mat3x3(cam_trf) * view_dir;
			ray.origin = cam_pos;

			const vec3 c_a(0.0f, 0.0f, 1.0f);
			const vec3 c_b(1.0f);
			vec3 color = mix(c_a, c_b, (f32)h / height);

			ray_hit hit;

			if(ray_cast_plane(ray, plane, 0.0f, 100.0f, &hit))
			{
				u32 sample_count = 5;
				u32 intersection_count = 0;
				for(u32 i = 0; i < sample_count; ++i)
				{
					vec3 rand_dir = rand_direction_in_cone(hit.normal, M_PI / 2.0f - 0.01f);

					struct ray new_ray;
					new_ray.direction = rand_dir;
					new_ray.origin = hit.point;
					struct ray_hit new_hit;
					bool intersects = ray_cast_sphere(new_ray, sphere, 0.01f, 100.0f, &new_hit);
					intersection_count += intersects;
				}

				color = vec3(1.0f, 0.0f, 0.0f) * (1.0f - ((f32)intersection_count / sample_count));
			}

			if(ray_cast(ray, sphere.base, 0.0f, 100.0f, &hit))
			{
				u32 sample_count = 5;
				u32 intersection_count = 0;
				for(u32 i = 0; i < sample_count; ++i)
				{
					vec3 rand_dir = rand_direction_in_cone(hit.normal, M_PI / 2.0f - 0.01f);

					struct ray new_ray;
					new_ray.direction = rand_dir;
					new_ray.origin = hit.point;
					struct ray_hit new_hit;
					bool intersects = ray_cast_plane(new_ray, plane, 0.01f, 100.0f, &new_hit);
					intersection_count += intersects;
				}

				color = vec3(0.0f, 1.0f, 0.0f) * (1.0f - ((f32)intersection_count / sample_count));
			}

			color = clamp(color, 0.0f, 1.0f);
			data[h][w][0] = color.r * 255;
			data[h][w][1] = color.g * 255;
			data[h][w][2] = color.b * 255;
		}

	save_tga(L"./image.tga", &data[0][0][0], sizeof(data), width, height);

	return 0;
}