#define _GNU_SOURCE 1
#include <engine/rt.h>
#include <pthread.h>

static void init_subsamples_8(u32_t width, u32_t height, vec2_t *subsamples) {
	const vec2_t SAMPLE_LOCS_8[8] = {vec2_init_2f(-7.0, 1.0), vec2_init_2f(-5.0, -5.0), vec2_init_2f(-1.0, -3.0),
	                                 vec2_init_2f(3.0, -7.0), vec2_init_2f(5.0, -1.0),  vec2_init_2f(7.0, 7.0),
	                                 vec2_init_2f(1.0, 3.0),  vec2_init_2f(-3.0, 5.0)};

	const vec2_t tex_size = 1.0f / vec2_init_2f(width, height); // Texel size

	for(u32_t i = 0; i < 8; ++i) {
		const vec2_t s = SAMPLE_LOCS_8[i] / 8.0f; // In [-1, 1]

		vec2_t sub_sample = s * tex_size; // In [-tex_size, tex_size]
		subsamples[i] = sub_sample;
	}
}

static void init_subsamples_n(u32_t width, u32_t height, vec2_t *subsamples, u32_t subsample_count) {
	const vec2_t tex_size = 1.0f / vec2_init_2f(width, height); // Texel size

	for(u32_t i = 0; i < subsample_count; ++i) {
		const vec2_t s = rand_0f_to_1f() * 2.0f - 1.0f; // In [-1, 1]
		vec2_t sub_sample = s * tex_size;               // In [-tex_size, tex_size]
		subsamples[i] = sub_sample;
	}
}

static mesh_t create_cube(vec3_t min, vec3_t max) {
	quadrilateral_t *quads = malloc(sizeof(quadrilateral_t) * 6);
	u32_t i = 0;
	quads[i++] = quad_init_yz(min.x, vec2_init_2f(min.y, max.y), vec2_init_2f(min.z, max.z), TRUE);  // Left
	quads[i++] = quad_init_yz(max.x, vec2_init_2f(min.y, max.y), vec2_init_2f(min.z, max.z), FALSE); // Right
	quads[i++] = quad_init_xz(vec2_init_2f(min.x, max.x), min.y, vec2_init_2f(min.z, max.z), TRUE);  // Bottom
	quads[i++] = quad_init_xz(vec2_init_2f(min.x, max.x), max.y, vec2_init_2f(min.z, max.z), FALSE); // Top
	quads[i++] = quad_init_xy(vec2_init_2f(min.x, max.x), vec2_init_2f(min.y, max.y), max.z, FALSE); // Front
	quads[i++] = quad_init_xy(vec2_init_2f(min.x, max.x), vec2_init_2f(min.y, max.y), min.z, TRUE);  // Back

	mesh_t mesh;
	memset(&mesh, 0, sizeof(mesh));
	mesh.quad_count = i;
	mesh.quads = quads;
	return mesh;
}

static vec3_t trace(const render_queue_t *rgraph, const ray_t *ray, u32_t depth, u32_t max_depth) {
	ray_hit_t hit;
	const material_t *hit_mtl;
	vec3_t color = {0};

	if(render_queue_closest_hit(rgraph, ray, &hit, &hit_mtl)) {
		vec3_t attenuation;
		pdf_t mtl_pdf;
		const vec3_t emitted = hit_mtl->emit_callback(hit_mtl, &hit);
		if(depth < max_depth && hit_mtl->scatter_callback(hit_mtl, ray, &hit, &attenuation, &mtl_pdf)) {
			pdf_t hittable_pdf = pdf_init_hittable(&rgraph->important_areas[0].box, hit.point);
			pdf_t mix_pdf = pdf_init_mixture(&hittable_pdf, &mtl_pdf);

			ray_t scatter_ray;
			scatter_ray.origin = hit.point;
			scatter_ray.direction = pdf_generate(&mix_pdf);
			const f32_t pdf_val = pdf_compute_value(&mix_pdf, scatter_ray.direction);

			color = attenuation * hit_mtl->scattering_pdf_callback(hit_mtl, ray, &hit, &scatter_ray);
			color *= trace(rgraph, &scatter_ray, depth + 1, max_depth) / pdf_val;
			color += emitted;
			color = vec3_de_nan(color);
		} else {
			color = emitted;
		}
	} else {
#if 0
		const vec3_t c_a = {0.9f, 0.9f, 1.0f};
		const vec3_t c_b = {0.4f, 0.6f, 1.0f};
		const f32_t t = 0.5f * (ray->direction.y + 1.0f);
		color = vec3_mix(c_a, c_b, t);
#else
		color = vec3_init_f(0.0f);
#endif
	}

	return color;
}

#if 0
static void init_sphere_lambertian(vec3_t pos, vec3_t prev_pos, f32_t radius, vec3_t color, renderable_t *r) {
	r->shape.sphere = sphere_init(vec3_init_f(0.0f), radius);
	r->shape_type = RENDERABLE_SHAPE_TYPE_SPHERE;
	r->world_transform = transform_init_t(pos);
	r->previous_world_transform = transform_init_t(prev_pos);
	r->material = material_init_lambertian();
	r->material.albedo_texture = texture_init_constant(color);
}

static void init_sphere_metal(vec3_t pos, vec3_t prev_pos, f32_t radius, vec3_t color, f32_t fuzz, renderable_t *r) {
	r->shape.sphere = sphere_init(vec3_init_f(0.0f), radius);
	r->shape_type = RENDERABLE_SHAPE_TYPE_SPHERE;
	r->world_transform = transform_init_t(pos);
	r->previous_world_transform = transform_init_t(prev_pos);
	r->material = material_init_metal();
	r->material.albedo_texture = texture_init_constant(color);
	r->material.metal_fuzz_texture = texture_init_constant(vec3_init_f(fuzz));
}

static void init_sphere_dielectric(vec3_t pos, vec3_t prev_pos, f32_t radius, f32_t refl_idx, renderable_t *r) {
	r->shape.sphere = sphere_init(vec3_init_f(0.0f), radius);
	r->shape_type = RENDERABLE_SHAPE_TYPE_SPHERE;
	r->world_transform = transform_init_t(pos);
	r->previous_world_transform = transform_init_t(prev_pos);
	r->material = material_init_dielectric();
	r->material.dielectric_reflection_index = texture_init_constant(vec3_init_f(refl_idx));
}

static render_queue_t random_scene() {
	const u32_t n = 500;
	renderable_t *spheres = (renderable_t *)malloc(sizeof(renderable_t) * (n + 1));
	memset(spheres, 0, sizeof(renderable_t) * (n + 1));

	init_sphere_lambertian(vec3_init_3f(0.0f, -1000.0f, 0.0f), vec3_init_3f(0.0f, -1000.0f, 0.0f), 1000.0f,
	                       vec3_init_f(0.0f), &spheres[0]);
	spheres[0].material.albedo_texture = texture_init_checker(vec3_init_3f(0.2f, 0.3f, 0.1f), vec3_init_f(0.9f));

	u32_t i = 1;
	for(i32_t a = -11; a < 11; a++) {
		for(i32_t b = -11; b < 11; b++) {
			const f32_t choose_mat = rand_0f_to_1f();
			const vec3_t center = vec3_init_3f(a + 0.9f * rand_0f_to_1f(), 0.2f, b + 0.9f * rand_0f_to_1f());
			if(vec3_length(center - vec3_init_3f(4.0f, 0.2f, 0.0f)) > 0.9f) {
				renderable_t *sphere = &spheres[i++];
				sphere->shape.sphere = sphere_init(center, 0.2f);

				const vec3_t rand_dir = random_in_unit_sphere();
				const f32_t rand_disp = rand_0f_to_1f() * 0.3f;
				const vec3_t prev_center = center + rand_dir * rand_disp;

				if(choose_mat < 0.8f) {
					const vec3_t color =
					    vec3_init_3f(rand_0f_to_1f() * rand_0f_to_1f(), rand_0f_to_1f() * rand_0f_to_1f(),
					                 rand_0f_to_1f() * rand_0f_to_1f());
					init_sphere_lambertian(center, prev_center, 0.2f, color, sphere);
				} else if(choose_mat < 0.95f) {
					const vec3_t color = vec3_init_3f(0.5f * (1.0f + rand_0f_to_1f()), 0.5f * (1.0f + rand_0f_to_1f()),
					                                  0.5f * (1.0f + rand_0f_to_1f()));

					init_sphere_metal(center, prev_center, 0.2f, color, 0.5f * rand_0f_to_1f(), sphere);
				} else {
					init_sphere_dielectric(center, prev_center, 0.2f, 1.5f, sphere);
				}
			}
		}
	}

	renderable_t *sphere = &spheres[i++];
	init_sphere_dielectric(vec3_init_3f(0.0f, 1.0f, 0.0f), vec3_init_3f(0.0f, 1.0f, 0.0f), 1.0f, 1.5f, sphere);

	sphere = &spheres[i++];
	init_sphere_lambertian(vec3_init_3f(-4.0f, 1.0f, 0.0f), vec3_init_3f(-4.0f, 1.0f, 0.0f), 1.0f, vec3_init_f(0.0f),
	                       sphere);
	sphere->material.albedo_texture = texture_init_image("data/earth.jpg");

	sphere = &spheres[i++];
	init_sphere_metal(vec3_init_3f(4.0f, 1.0f, 0.0f), vec3_init_3f(4.0f, 1.0f, 0.0f), 1.0f,
	                  vec3_init_3f(0.7f, 0.6f, 0.5f), 0.0f, sphere);

	assert(i <= n + 1);
	render_queue_t rgraph;
	rgraph.renderable_count = i;
	rgraph.renderables = spheres;
	return rgraph;
}
#endif

render_queue_t cornell_box() {
	const u32_t renderable_count = 8;
	renderable_t *renderables = malloc(sizeof(renderable_t) * renderable_count);
	memset(renderables, 0, sizeof(renderable_t) * renderable_count);
	u32_t i = 0;
	important_area_t *important_areas = malloc(sizeof(important_area_t) * 1);
	memset(important_areas, 0, sizeof(*important_areas) * 1);
	u32_t important_area_count = 0;

	renderable_t *light = &renderables[i++];
	light->material = material_init_emissive();
	light->material.emissive_texture = texture_init_constant(vec3_init_f(15.0f));
	light->shape_type = RENDERABLE_SHAPE_TYPE_MESH;
	light->shape.mesh.quad_count = 1;
	light->shape.mesh.quads = malloc(sizeof(quadrilateral_t) * 1);
	light->shape.mesh.quads[0] = quad_init_xz(vec2_init_2f(213.0f, 343.0f), 554.0f, vec2_init_2f(227.0f, 332.0f), TRUE);
	light->world_transform = light->previous_world_transform = transform_init_identity();

	important_area_t *light_important_area = &important_areas[important_area_count++];
	light_important_area->box =
	    aabb_init(vec3_init_3f(213.0f, 554.0f, 227.0f), vec3_init_3f(343.0f, 554.0f + EPSILON, 332.0f));

	renderable_t *left_wall = &renderables[i++];
	left_wall->material = material_init_lambertian();
	left_wall->material.albedo_texture = texture_init_constant(vec3_init_3f(0.12, 0.45, 0.15));
	left_wall->shape_type = RENDERABLE_SHAPE_TYPE_MESH;
	left_wall->shape.mesh.quad_count = 1;
	left_wall->shape.mesh.quads = malloc(sizeof(quadrilateral_t) * 1);
	left_wall->shape.mesh.quads[0] = quad_init_yz(555.0f, vec2_init_2f(0.0f, 555.0f), vec2_init_2f(0.0f, 555.0f), TRUE);
	left_wall->world_transform = left_wall->previous_world_transform = transform_init_identity();

	renderable_t *right_wall = &renderables[i++];
	right_wall->material = material_init_lambertian();
	right_wall->material.albedo_texture = texture_init_constant(vec3_init_3f(0.65, 0.05, 0.05));
	right_wall->shape_type = RENDERABLE_SHAPE_TYPE_MESH;
	right_wall->shape.mesh.quad_count = 1;
	right_wall->shape.mesh.quads = malloc(sizeof(quadrilateral_t) * 1);
	right_wall->shape.mesh.quads[0] = quad_init_yz(0.0f, vec2_init_2f(0.0f, 555.0f), vec2_init_2f(0.0f, 555.0f), FALSE);
	right_wall->world_transform = right_wall->previous_world_transform = transform_init_identity();

	renderable_t *top_wall = &renderables[i++];
	top_wall->material = material_init_lambertian();
	top_wall->material.albedo_texture = texture_init_constant(vec3_init_f(0.73));
	top_wall->shape_type = RENDERABLE_SHAPE_TYPE_MESH;
	top_wall->shape.mesh.quad_count = 1;
	top_wall->shape.mesh.quads = malloc(sizeof(quadrilateral_t) * 1);
	top_wall->shape.mesh.quads[0] = quad_init_xz(vec2_init_2f(0.0f, 555.0f), 555.0, vec2_init_2f(0.0f, 555.0f), TRUE);
	top_wall->world_transform = top_wall->previous_world_transform = transform_init_identity();

	renderable_t *bottom_wall = &renderables[i++];
	bottom_wall->material = material_init_lambertian();
	bottom_wall->material.albedo_texture = texture_init_constant(vec3_init_f(0.73));
	bottom_wall->shape_type = RENDERABLE_SHAPE_TYPE_MESH;
	bottom_wall->shape.mesh.quad_count = 1;
	bottom_wall->shape.mesh.quads = malloc(sizeof(quadrilateral_t) * 1);
	bottom_wall->shape.mesh.quads[0] = quad_init_xz(vec2_init_2f(0.0f, 555.0f), 0.0, vec2_init_2f(0.0f, 555.0f), FALSE);
	bottom_wall->world_transform = bottom_wall->previous_world_transform = transform_init_identity();

	renderable_t *back_wall = &renderables[i++];
	back_wall->material = material_init_lambertian();
	back_wall->material.albedo_texture = texture_init_constant(vec3_init_f(0.73));
	back_wall->shape_type = RENDERABLE_SHAPE_TYPE_MESH;
	back_wall->shape.mesh.quad_count = 1;
	back_wall->shape.mesh.quads = malloc(sizeof(quadrilateral_t) * 1);
	back_wall->shape.mesh.quads[0] = quad_init_xy(vec2_init_2f(0.0f, 555.0f), vec2_init_2f(0.0f, 555.0f), 555.0f, TRUE);
	back_wall->world_transform = back_wall->previous_world_transform = transform_init_identity();

	renderable_t *box1 = &renderables[i++];
	box1->material = material_init_lambertian();
	box1->material.albedo_texture = texture_init_constant(vec3_init_f(0.75f));
	box1->shape_type = RENDERABLE_SHAPE_TYPE_MESH;
	vec3_t min = vec3_init_3f(130, 0, 65);
	vec3_t max = vec3_init_3f(295, 160, 230);
	box1->shape.mesh = create_cube(-(max - min) / 2.0f, (max - min) / 2.0f);
	box1->world_transform = box1->previous_world_transform = transform_init_identity();
	box1->world_transform.translation = (min + max) / 2.0f;
	box1->world_transform.rotation = mat3_init_axis_angles(vec3_init_3f(0.0f, 1.0f, 0.0f), to_rad(-18.0f));

	renderable_t *tall_box = &renderables[i++];
	tall_box->material = material_init_lambertian();
	tall_box->material.albedo_texture = texture_init_constant(vec3_init_f(0.75));
	tall_box->shape_type = RENDERABLE_SHAPE_TYPE_MESH;
	min = vec3_init_3f(265, 0, 295);
	max = vec3_init_3f(430, 330, 460);
	tall_box->shape.mesh = create_cube(-(max - min) / 2.0f, (max - min) / 2.0f);
	tall_box->world_transform = tall_box->previous_world_transform = transform_init_identity();
	tall_box->world_transform.translation = (min + max) / 2.0f;
	tall_box->world_transform.rotation = mat3_init_axis_angles(vec3_init_3f(0.0f, 1.0f, 0.0f), to_rad(15.0f));

	render_queue_t rgraph = {
	    0,
	};
	rgraph.renderable_count = i;
	rgraph.renderables = renderables;
	rgraph.important_areas = important_areas;
	rgraph.important_area_count = important_area_count;
	return rgraph;
}

typedef struct run_context_t {
	u32_t tile_size;
	u32_t width;
	u32_t height;
	u32_t max_trace_depth;
	vec2_t *subsamples;
	u32_t subsample_count;
	mat4_t camera_transform;
	mat4_t invert_vp_matrix;
	vec3_t camera_position;
	render_queue_t render_graph;
	u8_t *pixel_buffer;
	u32_t last_tile;
} run_context_t;

static void *run_thread(void *user_data) {
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

				for(u32_t s = 0; s < ctx->subsample_count; ++s) {
					const vec2_t s_ndc = ndc + ctx->subsamples[s];

					const vec4_t view_dir4 =
					    mat4_mul_vec4(&ctx->invert_vp_matrix, vec4_init_4f(s_ndc.x, s_ndc.y, 1.0f, 1.0f));
					const vec3_t far_point = view_dir4.xyz / view_dir4.w;

					const ray_t primary_ray =
					    ray_init(ctx->camera_position, vec3_normalize(far_point - ctx->camera_position));

					color += trace(&ctx->render_graph, &primary_ray, 0, ctx->max_trace_depth);
				}

				color /= (f32_t)ctx->subsample_count;

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

	const u32_t width = 500;
	const u32_t height = 500;
	const u32_t subsample_count = 16;

	run_context_t ctx;
	memset(&ctx, 0, sizeof(ctx));
	ctx.tile_size = 64;
	ctx.width = width;
	ctx.height = height;
	ctx.max_trace_depth = 2;
	vec2_t subsamples[subsample_count];
	ctx.subsamples = subsamples;
	ctx.subsample_count = subsample_count;
	if(subsample_count == 8) {
		init_subsamples_8(width, height, ctx.subsamples);
	} else {
		init_subsamples_n(width, height, ctx.subsamples, subsample_count);
	}

	// const vec3_t cam_pos = vec3_init_3f(0.0f, 1.1f, 6.0f);
	const vec3_t cam_pos = vec3_init_3f(278.0f, 278.0f, -800.0f);
	// const vec3_t ref_point = cam_pos + vec3_init_3f(0.0f, -0.1f, -1.0f);
	const vec3_t ref_point = cam_pos + vec3_init_3f(0.0f, 0.0f, 1.0f);
	// const vec3_t ref_point = vec3_init_3f(278.0f, 278.0f, 0.0f);
	const vec3_t up = vec3_init_3f(0.0f, 1.0f, 0.0f);
	const mat4_t view_mat = look_at(cam_pos, ref_point, up);
	const mat4_t cam_trf = mat4_invert(&view_mat);
	const mat4_t proj_mat = perspective(to_rad(40.0f), (f32_t)width / (f32_t)height, 0.1f, 1000.0f);
	const mat4_t vp_mat = mat4_mul_mat4(&proj_mat, &view_mat);
	const mat4_t inv_vp_mat = mat4_invert(&vp_mat);

	ctx.camera_transform = cam_trf;
	ctx.invert_vp_matrix = inv_vp_mat;
	ctx.camera_position = cam_pos;

	// Render graph
#if 0
	renderable_t s[5];
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
	// ctx.render_graph = random_scene();
	ctx.render_graph = cornell_box();
#endif

	u8_t *pixel_buffer = malloc(height * width * 3);
	ctx.pixel_buffer = pixel_buffer;

	const u32_t thread_count = get_gpu_core_count();
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

	save_tga("./image.tga", pixel_buffer, height * width * 3, width, height);

	return 0;
}
