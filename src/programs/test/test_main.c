#include <engine/rt.h>

void init_subsamples(u32_t width, u32_t height, vec2_t *subsamples) {
	vec2_t SAMPLE_LOCS_8[8] = {vec2_init_2f(-7.0, 1.0), vec2_init_2f(-5.0, -5.0), vec2_init_2f(-1.0, -3.0),
	                           vec2_init_2f(3.0, -7.0), vec2_init_2f(5.0, -1.0),  vec2_init_2f(7.0, 7.0),
	                           vec2_init_2f(1.0, 3.0),  vec2_init_2f(-3.0, 5.0)};

	vec2_t tex_size = 1.0f / vec2_init_2f(width, height); // Texel size

	for(u32_t i = 0; i < 8; ++i) {
		const vec2_t s = SAMPLE_LOCS_8[i] / 8.0f; // In [-1, 1]

		vec2_t sub_sample = s * tex_size; // In [-tex_size, tex_size]
		subsamples[i] = sub_sample;
	}
}

bool_t lambertian_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                          ray_t *scattered_ray) {
	(void)in_ray;
	const vec3_t target = hit->point + hit->normal + random_in_init_sphere();
	*scattered_ray = ray_init(hit->point, vec3_normalize(target - hit->point));
	*attenuation = mtl->albedo;
	return TRUE;
}

bool_t metal_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                     ray_t *scattered_ray) {
	const vec3_t reflected = vec3_reflect(in_ray->direction, hit->normal);
	*scattered_ray = ray_init(hit->point, vec3_normalize(reflected));
	*attenuation = mtl->albedo;
	return vec3_dot(reflected, hit->normal) > 0.0f;
}

bool_t closest_hit(const render_graph_t *rgraph, const ray_t *ray, ray_hit_t *closest_hit, const material_t **hit_mtl) {
	closest_hit->t = INFINITY;
	bool_t has_hit = FALSE;

	for(u32_t i = 0; i < rgraph->sphere_count; ++i) {
		ray_hit_t hit;
		if(ray_cast_sphere(ray, &rgraph->spheres[i].sphere, 0.001f, 1000.f, &hit)) {
			if(hit.t < closest_hit->t) {
				*closest_hit = hit;
				*hit_mtl = &rgraph->spheres[i].material;
				has_hit = TRUE;
			}
		}
	}

	return has_hit;
}

vec3_t trace(const render_graph_t *rgraph, const ray_t *ray, u32_t depth) {
	ray_hit_t hit;
	const material_t *hit_mtl;
	vec3_t color = {0};

	if(closest_hit(rgraph, ray, &hit, &hit_mtl)) {
		ray_t new_ray;
		vec3_t attenuation;
		if(depth < 50 && hit_mtl->scatter_callback(hit_mtl, ray, &hit, &attenuation, &new_ray)) {
			color = attenuation * trace(rgraph, &new_ray, depth + 1);
		} else {
			color = vec3_init_f(0.0f);
		}
	} else {
		const vec3_t c_a = {0.9f, 0.9f, 1.0f};
		const vec3_t c_b = {0.4f, 0.6f, 1.0f};
		const f32_t t = 0.5f * (ray->direction.y + 1.0f);
		color = vec3_mix(c_a, c_b, t);
	}

	return color;
}

int main(int argc, char **argv) {
	(void)argc;
	(void)argv;
	const u32_t width = 1024;
	const u32_t height = 768;

	const vec3_t cam_pos = {0.0f, 0.5f, 10.0f};
	const vec3_t eye = {0.0f, 0.0f, -1.0f};
	const vec3_t up = {0.0f, 1.0f, 0.0f};
	const mat4_t view_mat = look_at(eye, cam_pos, up);
	const mat4_t cam_trf = mat4_invert(&view_mat);
	const mat4_t proj_mat = perspective(to_rad(60.0f), (f32_t)width / (f32_t)height, 0.1f, 100.0f);
	const mat4_t vp_mat = mat4_mul_mat4(&proj_mat, &view_mat);
	const mat4_t inv_vp_mat = mat4_invert(&vp_mat);

	// Subpixels
	vec2_t subsamples[8];
	init_subsamples(width, height, subsamples);

	// Render graph
	renderable_sphere_t s[4] = {0};
	s[0].sphere = sphere_init(vec3_init_3f(0.0f, 0.0f, -1.0f), 0.5f);
	s[0].material.scatter_callback = lambertian_scatter;
	s[0].material.albedo = vec3_init_3f(0.8f, 0.3f, 0.3f);
	s[1].sphere = sphere_init(vec3_init_3f(0.0f, -100.5f, -1.0f), 100.0f);
	s[1].material.scatter_callback = lambertian_scatter;
	s[1].material.albedo = vec3_init_3f(0.8f, 0.8f, 0.0f);
	s[2].sphere = sphere_init(vec3_init_3f(1.0f, 0.0f, -1.0f), 0.5f);
	s[2].material.scatter_callback = metal_scatter;
	s[2].material.albedo = vec3_init_3f(0.8f, 0.6f, 0.2f);
	s[3].sphere = sphere_init(vec3_init_3f(-1.0f, 0.0f, -1.0f), 0.5f);
	s[3].material.scatter_callback = metal_scatter;
	s[3].material.albedo = vec3_init_3f(0.8f, 0.8f, 0.8f);
	render_graph_t rgraph;
	rgraph.spheres = s;
	rgraph.sphere_count = NELEMS(s);

	u8_t data[height][width][3];
	for(u32_t w = 0; w < width; ++w) {
		for(u32_t h = 0; h < height; ++h) {
			vec2_t ndc = {(w + 0.5f) / width, (h + 0.5f) / height};
			ndc = ndc * vec2_init_f(2.0f) - vec2_init_f(1.0f);

			vec3_t color = vec3_init_f(0.0f);

			for(u32_t s = 0; s < 8; ++s) {
				const vec2_t s_ndc = ndc + subsamples[s];

				const vec3_t view_dir = mat4_mul_vec4(&inv_vp_mat, vec4_init_4f(s_ndc.x, s_ndc.y, 1.0f, 1.0f)).xyz;
				mat3_t cam_rot = mat3_init_mat4(&cam_trf);
				ray_t primary_ray;
				primary_ray.direction = mat3_mul_vec3(&cam_rot, view_dir);
				primary_ray.origin = cam_pos;

				color += trace(&rgraph, &primary_ray, 0);
			}

			color /= vec3_init_f(NELEMS(subsamples));

			// Write
			const f32_t gamma = 2.2f;
			color = vec3_init_3f(powf(color.x, 1.0f / gamma), powf(color.y, 1.0f / gamma), powf(color.z, 1.0f / gamma));
			color = vec3_clamp(color, vec3_init_f(0.0f), vec3_init_f(1.0f));
			data[h][w][0] = color.x * 255;
			data[h][w][1] = color.y * 255;
			data[h][w][2] = color.z * 255;
		}
	}

	save_tga("./image.tga", &data[0][0][0], sizeof(data), width, height);

	return 0;
}
