#include <engine/misc/public/tga.hpp>
#include <engine/math/public/all.hpp>
#include <engine/util/public/all.hpp>
#include <engine/renderer/public/render_queue.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <cstdio>

int main(int, char**)
{
	using namespace RT_NAMESPACE;

	const u32_t width = 512;
	const u32_t height = 256;

	vec3_t cam_pos = vec3_t(0.0f, 1.0f, 0.0f);

	mat4_t view_mat = lookAt(vec3_t(0.0f, 0.0f, -1.0f), cam_pos, vec3_t(0.0f, 1.0f, 0.0f));
	mat4_t cam_trf = inverse(view_mat);
	mat4_t proj_mat = perspective(glm::radians(45.0f), (f32_t)width / (f32_t)height, 0.1f, 100.0f);
	mat4_t vp_mat = proj_mat * view_mat;
	mat4_t inv_vp_mat = inverse(vp_mat);

	render_queue_t rqueue;
	rqueue.renderables = new_array<renderable_t>(2);
	rqueue.renderable_count = 2;

	rqueue.renderables[0].sphere.center = vec3_t(0.0f, 1.0f, -5.0f);
	rqueue.renderables[0].sphere.radius = 1.0f;

	rqueue.renderables[1].plane.offset = 0.0f;
	rqueue.renderables[1].plane.normal = vec3_t(0.0f, 1.0f, 0.0f);

	u8_t data[height][width][3];
	for(u32_t w = 0; w < width; ++w)
		for(u32_t h = 0; h < height; ++h)
		{
			vec2_t ndc((w + 0.5f) / width, (h + 0.5f) / height);
			ndc = ndc * 2.0f - 1.0f;

			vec3_t view_dir = inv_vp_mat * vec4_t(ndc, 1.0f, 1.0f);
			ray_t ray;
			ray.direction = mat3_t(cam_trf) * view_dir;
			ray.origin = cam_pos;

			const vec3_t c_a(0.0f, 0.0f, 1.0f);
			const vec3_t c_b(1.0f);
			vec3_t color = mix(c_a, c_b, (f32_t)h / height);

			ray_hit_t hit;

			/*if(ray_cast_plane(ray, plane, 0.0f, 100.0f, hit))
			{
				u32 sample_count = 5;
				u32 intersection_count = 0;
				for(u32 i = 0; i < sample_count; ++i)
				{
					vec3_t rand_dir = rand_direction_in_cone(hit.normal, M_PI / 2.0f - 0.01f);

					struct ray new_ray;
					new_ray.direction = rand_dir;
					new_ray.origin = hit.point;
					struct ray_hit new_hit;
					boolean intersects = ray_cast_sphere(new_ray, sphere, 0.01f, 100.0f, new_hit);
					intersection_count += (u32)intersects;
				}

				color = vec3_t(1.0f, 0.0f, 0.0f) * (1.0f - ((f32)intersection_count / sample_count));
			}

			if(ray_cast(ray, sphere.base, 0.0f, 100.0f, hit))
			{
				u32 sample_count = 5;
				u32 intersection_count = 0;
				for(u32 i = 0; i < sample_count; ++i)
				{
					vec3_t rand_dir = rand_direction_in_cone(hit.normal, M_PI / 2.0f - 0.01f);

					struct ray new_ray;
					new_ray.direction = rand_dir;
					new_ray.origin = hit.point;
					struct ray_hit new_hit;
					boolean intersects = ray_cast_plane(new_ray, plane, 0.01f, 100.0f, new_hit);
					intersection_count += (u32)intersects;
				}

				color = vec3_t(0.0f, 1.0f, 0.0f) * (1.0f - ((f32)intersection_count / sample_count));
			}*/

			color = clamp(color, 0.0f, 1.0f);
			data[h][w][0] = color.r * 255;
			data[h][w][1] = color.g * 255;
			data[h][w][2] = color.b * 255;
		}

	save_tga(L"./image.tga", &data[0][0][0], sizeof(data), width, height);

	return 0;
}
