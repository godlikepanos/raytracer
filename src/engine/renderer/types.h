#pragma once

#include <engine/math.h>

struct material_t;
struct texture_t;

typedef bool_t (*scatter_callback_t)(const struct material_t *mtl, const ray_t *ray, const ray_hit_t *hit,
                                     vec3_t *attenuation, ray_t *scattered_ray);

typedef vec3_t (*texture_callback_t)(const struct texture_t *tex, vec2_t uv, vec3_t point);

typedef struct texture_t {
	vec3_t color0;
	vec3_t color1;
	texture_callback_t callback;
} texture_t;

typedef struct material_t {
	scatter_callback_t scatter_callback;
	texture_t albedo_texture;
	texture_t metal_fuzz_texture;
	texture_t dielectric_reflection_index;
} material_t;

typedef struct renderable_t {
	union {
		sphere_t sphere;
		aabb_t box;
	} shape;
	vec3_t previous_position;
	material_t material;
} renderable_t;

typedef struct render_queue_t {
	const renderable_t *spheres;
	u32_t sphere_count;
} render_queue_t;

typedef struct bvh_node_t {
	aabb_t box;
	struct bvh_node_t *left;
	struct bvh_node_t *right;
	u32_t *sphere_list;
	u32_t sphere_count;
} bvh_node_t;

typedef struct bvh_tree_t {
	bvh_node_t root;
} bvh_tree_t;
