#pragma once

#include <engine/util/std.h>
#include <math.h>

static const f32_t PI = 3.14159265358979323846264338327950288f;
static const f32_t EPSILON = 1.0e-6f;

typedef float vec2_t __attribute__((ext_vector_type(2)));
typedef float vec3_t __attribute__((ext_vector_type(3)));
typedef float vec4_t __attribute__((ext_vector_type(4)));

typedef struct mat3_t {
	vec3_t m[3];
} mat3_t;

typedef struct mat4_t {
	vec4_t m[4];
} mat4_t;

typedef enum collision_shape_type_e {
	COLLISION_SHAPE_TYPE_SPHERE,
	COLLISION_SHAPE_TYPE_PLANE,
} collision_shape_type_e;

typedef struct ray_t {
	vec3_t origin;
	vec3_t direction;
} ray_t;

typedef struct ray_hit_t {
	vec3_t point;
	vec3_t normal;
	f32_t t;
	vec2_t uv;
} ray_hit_t;

typedef struct sphere_t {
	vec3_t center;
	f32_t radius;
} sphere_t;

typedef struct plane_t {
	vec3_t normal;
	f32_t offset;
} plane_t;

typedef struct aabb_t {
	vec3_t min;
	vec3_t max;
} aabb_t;

typedef struct triangle_t {
	vec3_t vertices[3];
} triangle_t;

typedef struct quadrilateral_t {
	vec3_t vertices[4];
} quadrilateral_t;
