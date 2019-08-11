#pragma once

#include <engine/math.h>
#include <engine/renderer/types.h>
#include <external/stb_image.h>

// texture_t
static inline vec3_t texture_constant(const struct texture_t *tex, vec2_t uv, vec3_t point) {
	(void)uv;
	(void)point;
	return tex->color0;
}

static inline vec3_t texture_checker(const struct texture_t *tex, vec2_t uv, vec3_t point) {
	(void)uv;
	const f32_t sines = sin(10.0f * point.x) * sin(10.0f * point.y) * sin(10.0f * point.z);
	return (sines < 0.0f) ? tex->color0 : tex->color1;
}

static inline vec3_t texture_image(const struct texture_t *tex, vec2_t uv, vec3_t point) {
	(void)point;
	const u32_t x = CLAMP(uv.x * tex->image_width, 0.0f, (f32_t)tex->image_width - 1.0f);
	const u32_t y = CLAMP((1.0f - uv.y) * tex->image_height, 0.0f, (f32_t)tex->image_height - 1.0f);
	assert(x < tex->image_width && y < tex->image_height);
	const u32_t comp_count = 4;
	const f32_t r = tex->image_pixels[y * comp_count * tex->image_width + x * comp_count + 0];
	const f32_t g = tex->image_pixels[y * comp_count * tex->image_width + x * comp_count + 1];
	const f32_t b = tex->image_pixels[y * comp_count * tex->image_width + x * comp_count + 2];
	return vec3_init_3f(r, g, b) / 255.0f;
}

static inline texture_t texture_init_constant(vec3_t color) {
	texture_t tex;
	tex.color0 = color;
	tex.color1 = color;
	tex.callback = texture_constant;
	return tex;
}

static inline texture_t texture_init_checker(vec3_t color0, vec3_t color1) {
	texture_t tex;
	tex.color0 = color0;
	tex.color1 = color1;
	tex.callback = texture_checker;
	return tex;
}

static inline texture_t texture_init_image(const char_t *filename) {
	texture_t tex;
	int nx, ny, nn;
	u8_t *tex_data = stbi_load(filename, &nx, &ny, &nn, 4);
	if(tex_data == NULL) {
		fprintf(stderr, "Failed to open: %s\n", filename);
		tex = texture_init_constant(vec3_init_f(0.0f));
	} else {
		tex.callback = texture_image;
		tex.image_pixels = tex_data;
		tex.image_width = (u32_t)nx;
		tex.image_height = (u32_t)ny;
	}
	return tex;
}

// material_t
bool_t lambertian_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                          ray_t *scattered_ray);

bool_t metal_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                     ray_t *scattered_ray);

bool_t dielectric_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                          ray_t *scattered_ray);

static inline material_t material_init_lambertian() {
	material_t mtl;
	memset(&mtl, 0, sizeof(mtl));
	mtl.scatter_callback = lambertian_scatter;
	return mtl;
}

static inline material_t material_init_metal() {
	material_t mtl;
	memset(&mtl, 0, sizeof(mtl));
	mtl.scatter_callback = metal_scatter;
	return mtl;
}

static inline material_t material_init_dielectric() {
	material_t mtl;
	memset(&mtl, 0, sizeof(mtl));
	mtl.scatter_callback = dielectric_scatter;
	return mtl;
}

// render_queue_t
aabb_t render_queue_compute_aabb(const render_queue_t *rqueue);
