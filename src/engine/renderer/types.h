#pragma once

#include <engine/math.h>

struct material_t;

typedef bool_t (*scatter_callback_t)(const struct material_t *mtl, const ray_t *ray, const ray_hit_t *hit,
                                     vec3_t *attenuation, ray_t *scattered_ray);

typedef struct material_t {
	scatter_callback_t scatter_callback;
	vec3_t albedo;
	vec3_t metal_fuzz;
	f32_t dielectric_refl_idx;
} material_t;

typedef struct renderable_sphere_t {
	sphere_t sphere;
	material_t material;
} renderable_sphere_t;

typedef struct render_graph_t {
	const renderable_sphere_t *spheres;
	u32_t sphere_count;
} render_graph_t;
