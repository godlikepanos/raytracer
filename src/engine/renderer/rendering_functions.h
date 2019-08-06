#pragma once

#include <engine/math.h>
#include <engine/renderer/types.h>

static inline f32_t schlick(f32_t cosine, f32_t ref_idx) {
	f32_t r0 = (1.0f - ref_idx) / (1.0f + ref_idx);
	r0 *= r0;
	return r0 + (1.0f - r0) * powf(1.0f - cosine, 5.0f);
}

bool_t lambertian_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                          ray_t *scattered_ray);

bool_t metal_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                     ray_t *scattered_ray);

bool_t dielectric_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                          ray_t *scattered_ray);
