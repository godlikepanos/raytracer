#pragma once

#include <engine/math.h>
#include <engine/renderer/types.h>

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
