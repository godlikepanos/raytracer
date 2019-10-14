#pragma once

#include <engine/math.h>

struct material_t;
struct texture_t;
struct pdf_t;

// ** PDF **
typedef enum {
	PDF_TYPE_COSINE,
	PDF_TYPE_HITTABLE,
	PDF_TYPE_MIXTURE,
} pdf_type_e;

typedef struct pdf_t {
	union {
		// Cosine PDF
		struct {
			mat3_t orthonormal_basis;
		} cosine;

		// Hittable PDF
		struct {
			aabb_t box;
			vec3_t ray_origin;
		} hittable;

		// Mixture PDF
		struct {
			const struct pdf_t *pdfs[8];
			u32_t pdf_count;
		} mixture;
	};

	pdf_type_e pdf_type;
} pdf_t;

// ** Textures **
typedef vec3_t (*texture_callback_t)(const struct texture_t *tex, vec2_t uv, vec3_t point);

typedef struct texture_t {
	vec3_t color0;
	vec3_t color1;
	u8_t *image_pixels;
	u32_t image_width;
	u32_t image_height;
	texture_callback_t callback;
} texture_t;

// ** Materials **
typedef bool_t (*scatter_callback_t)(const struct material_t *mtl, const ray_t *ray, const ray_hit_t *hit,
                                     vec3_t *attenuation, ray_t *scattered_ray, f32_t *pdf);

typedef f32_t (*scatter_pdf_callback_t)(const struct material_t *mtl, const ray_t *ray, const ray_hit_t *hit,
                                        const ray_t *scattered_ray);

typedef vec3_t (*emit_callback_t)(const struct material_t *mtl, const ray_hit_t *hit);

typedef struct material_t {
	scatter_callback_t scatter_callback;
	emit_callback_t emit_callback;
	scatter_pdf_callback_t scatter_pdf_callback;
	texture_t albedo_texture;
	texture_t metal_fuzz_texture;
	texture_t dielectric_reflection_index;
	texture_t emissive_texture;
} material_t;

// ** Render Queue **
typedef struct mesh_t {
	triangle_t *triangles;
	quadrilateral_t *quads;
	u32_t triangle_count;
	u32_t quad_count;
} mesh_t;

typedef enum renderable_shape_type_e {
	RENDERABLE_SHAPE_TYPE_SPHERE,
	RENDERABLE_SHAPE_TYPE_AABB,
	RENDERABLE_SHAPE_TYPE_MESH,
	RENDERABLE_SHAPE_TYPE_COUNT
} renderable_shape_type_e;

typedef struct renderable_t {
	union {
		sphere_t sphere;
		aabb_t box;
		mesh_t mesh;
	} shape;
	renderable_shape_type_e shape_type;
	transform_t world_transform;
	transform_t previous_world_transform;
	material_t material;
} renderable_t;

typedef struct important_area_t {
	aabb_t box;
} important_area_t;

typedef struct render_queue_t {
	const renderable_t *renderables;
	u32_t renderable_count;

	const important_area_t *important_areas;
	u32_t important_area_count;
} render_queue_t;

// ** BVH **
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
