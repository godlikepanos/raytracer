#define _GNU_SOURCE 1
#include <engine/rt.h>
#include <pthread.h>

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

render_graph_t random_scene() {
	const u32_t n = 500;
	renderable_sphere_t *spheres = (renderable_sphere_t *)malloc(sizeof(renderable_sphere_t) * (n + 1));
	spheres[0].sphere = sphere_init(vec3_init_3f(0.0f, -1000.0f, 0.0f), 1000.0f);
	spheres[0].material.albedo = vec3_init_f(0.5f);
	spheres[0].material.scatter_callback = lambertian_scatter;
	u32_t i = 1;
	for(i32_t a = -11; a < 11; a++) {
		for(i32_t b = -11; b < 11; b++) {
			f32_t choose_mat = rand_0f_to_1f();
			vec3_t center = vec3_init_3f(a + 0.9f * rand_0f_to_1f(), 0.2f, b + 0.9f * rand_0f_to_1f());
			if(vec3_length(center - vec3_init_3f(4.0f, 0.2f, 0.0f)) > 0.9f) {
				renderable_sphere_t *sphere = &spheres[i++];
				sphere->sphere = sphere_init(center, 0.2f);
				if(choose_mat < 0.8f) {
					sphere->material.albedo =
					    vec3_init_3f(rand_0f_to_1f() * rand_0f_to_1f(), rand_0f_to_1f() * rand_0f_to_1f(),
					                 rand_0f_to_1f() * rand_0f_to_1f());
					sphere->material.scatter_callback = lambertian_scatter;
				} else if(choose_mat < 0.95f) {
					sphere->material.albedo =
					    vec3_init_3f(0.5f * (1.0f + rand_0f_to_1f()), 0.5f * (1.0f + rand_0f_to_1f()),
					                 0.5f * (1.0f + rand_0f_to_1f()));
					sphere->material.metal_fuzz = 0.5f * rand_0f_to_1f();
					sphere->material.scatter_callback = metal_scatter;
				} else {
					sphere->material.dielectric_refl_idx = 1.5f;
					sphere->material.scatter_callback = dielectric_scatter;
				}
			}
		}
	}

	renderable_sphere_t *sphere = &spheres[i++];
	sphere->sphere = sphere_init(vec3_init_3f(0.0f, 1.0f, 0.0f), 1.0f);
	sphere->material.dielectric_refl_idx = 1.5f;
	sphere->material.scatter_callback = dielectric_scatter;

	sphere = &spheres[i++];
	sphere->sphere = sphere_init(vec3_init_3f(-4.0f, 1.0f, 0.0f), 1.0f);
	sphere->material.albedo = vec3_init_3f(0.4f, 0.2f, 0.1f);
	sphere->material.scatter_callback = lambertian_scatter;

	sphere = &spheres[i++];
	sphere->sphere = sphere_init(vec3_init_3f(4.0f, 1.0f, 0.0f), 1.0f);
	sphere->material.albedo = vec3_init_3f(0.7f, 0.6f, 0.5f);
	sphere->material.scatter_callback = metal_scatter;
	sphere->material.metal_fuzz = 0.0f;
	assert(i <= n + 1);

	render_graph_t rgraph;
	rgraph.sphere_count = i;
	rgraph.spheres = spheres;
	return rgraph;
}

typedef struct run_context_t {
	u32_t tile_size;
	u32_t width;
	u32_t height;
	vec2_t subsamples[8];
	mat4_t camera_transform;
	mat4_t invert_vp_matrix;
	vec3_t camera_position;
	render_graph_t render_graph;
	u8_t *pixel_buffer;
	u32_t last_tile;
} run_context_t;

void *run_thread(void *user_data) {
	run_context_t *ctx = (run_context_t *)user_data;

	const u32_t tile_count_x = (ctx->width + ctx->tile_size - 1) / ctx->tile_size;
	const u32_t tile_count_y = (ctx->height + ctx->tile_size - 1) / ctx->tile_size;
	const u32_t tile_count = tile_count_x * tile_count_y;
	u32_t tile_idx;
	while((tile_idx = __atomic_fetch_add(&ctx->last_tile, 1, __ATOMIC_RELAXED)) < tile_count) {
		const u32_t tile_idx_x = tile_idx % tile_count_x;
		const u32_t tile_idx_y = tile_idx / tile_count_x;

		for(u32_t h = tile_idx_y * ctx->tile_size; h < MIN(ctx->height, (tile_idx_y + 1) * ctx->tile_size); ++h) {
			for(u32_t w = tile_idx_x * ctx->tile_size; w < MIN(ctx->width, (tile_idx_x + 1) * ctx->tile_size); ++w) {
				vec2_t ndc = {(w + 0.5f) / ctx->width, (h + 0.5f) / ctx->height};
				ndc = ndc * vec2_init_f(2.0f) - vec2_init_f(1.0f);

				vec3_t color = vec3_init_f(0.0f);

				for(u32_t s = 0; s < 8; ++s) {
					const vec2_t s_ndc = ndc + ctx->subsamples[s];

					const vec4_t view_dir4 =
					    mat4_mul_vec4(&ctx->invert_vp_matrix, vec4_init_4f(s_ndc.x, s_ndc.y, 1.0f, 1.0f));
					const vec3_t view_dir = view_dir4.xyz / view_dir4.w;

					mat3_t cam_rot = mat3_init_mat4(&ctx->camera_transform);
					const ray_t primary_ray =
					    ray_init(ctx->camera_position, vec3_normalize(mat3_mul_vec3(&cam_rot, view_dir)));

					color += trace(&ctx->render_graph, &primary_ray, 0);
				}

				color /= vec3_init_f(NELEMS(ctx->subsamples));

				// Write
				const f32_t gamma = 2.2f;
				color =
				    vec3_init_3f(powf(color.x, 1.0f / gamma), powf(color.y, 1.0f / gamma), powf(color.z, 1.0f / gamma));
				color = vec3_clamp(color, vec3_init_f(0.0f), vec3_init_f(1.0f));

				const u32_t pixel_buffer_idx = h * ctx->width * 3 + w * 3;
				assert(pixel_buffer_idx + 3 <= ctx->width * ctx->height * 3);
				ctx->pixel_buffer[pixel_buffer_idx + 0] = color.x * 255;
				ctx->pixel_buffer[pixel_buffer_idx + 1] = color.y * 255;
				ctx->pixel_buffer[pixel_buffer_idx + 2] = color.z * 255;
			}
		}
	}

	return NULL;
}

int main(int argc, char **argv) {
	(void)argc;
	(void)argv;

	seed_mt(time(NULL));
	srand(time(NULL));

	const u32_t width = 1920;
	const u32_t height = 1080;

	run_context_t ctx;
	memset(&ctx, 0, sizeof(ctx));
	ctx.tile_size = 64;
	ctx.width = width;
	ctx.height = height;
	init_subsamples(width, height, ctx.subsamples);

	const vec3_t cam_pos = vec3_init_3f(0.0f, 1.1f, 6.0f);
	const vec3_t ref_point = cam_pos + vec3_init_3f(0.0f, -0.1f, -1.0f);
	const vec3_t up = vec3_init_3f(0.0f, 1.0f, 0.0f);
	const mat4_t view_mat = look_at(cam_pos, ref_point, up);
	const mat4_t cam_trf = mat4_invert(&view_mat);
	const mat4_t proj_mat = perspective(to_rad(70.0f), (f32_t)width / (f32_t)height, 0.1f, 1000.0f);
	const mat4_t vp_mat = mat4_mul_mat4(&proj_mat, &view_mat);
	const mat4_t inv_vp_mat = mat4_invert(&vp_mat);

	ctx.camera_transform = cam_trf;
	ctx.invert_vp_matrix = inv_vp_mat;
	ctx.camera_position = cam_pos;

	// Render graph
#if 0
	renderable_sphere_t s[5];
	memset(s, 0, sizeof(s));
	s[0].sphere = sphere_init(vec3_init_3f(0.0f, 0.0f, -1.0f), 0.5f);
	s[0].material.scatter_callback = lambertian_scatter;
	s[0].material.albedo = vec3_init_3f(0.1f, 0.2f, 0.5f);
	s[1].sphere = sphere_init(vec3_init_3f(0.0f, -100.5f, -1.0f), 100.0f);
	s[1].material.scatter_callback = lambertian_scatter;
	s[1].material.albedo = vec3_init_3f(0.8f, 0.8f, 0.0f);
	s[2].sphere = sphere_init(vec3_init_3f(1.0f, 0.0f, -1.0f), 0.5f);
	s[2].material.scatter_callback = metal_scatter;
	s[2].material.albedo = vec3_init_3f(0.8f, 0.6f, 0.2f);
	s[2].material.metal_fuzz = 0.05f;
	s[3].sphere = sphere_init(vec3_init_3f(-1.0f, 0.0f, -1.0f), 0.5f);
	s[3].material.scatter_callback = dielectric_scatter;
	s[3].material.dielectric_refl_idx = 1.5f;
	s[4].sphere = sphere_init(vec3_init_3f(-1.0f, 0.0f, -1.0f), -0.45f);
	s[4].material.scatter_callback = dielectric_scatter;
	s[4].material.dielectric_refl_idx = 1.5f;
	ctx.render_graph.spheres = s;
	ctx.render_graph.sphere_count = NELEMS(s);
#else
	ctx.render_graph = random_scene();
#endif

	u8_t pixel_buffer[height * width * 3];
	ctx.pixel_buffer = pixel_buffer;

	const u32_t thread_count = 32;
	pthread_t threads[thread_count];
	for(u32_t i = 0; i < thread_count; ++i) {
		pthread_attr_t attr;
		cpu_set_t cpus;
		pthread_attr_init(&attr);

		CPU_ZERO(&cpus);
		CPU_SET(i, &cpus);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);

		pthread_create(&threads[i], &attr, run_thread, &ctx);
	}

	for(u32_t i = 0; i < thread_count; ++i) {
		void *out;
		pthread_join(threads[i], &out);
	}

	save_tga("./image.tga", pixel_buffer, sizeof(pixel_buffer), width, height);

	return 0;
}
