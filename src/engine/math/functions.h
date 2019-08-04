#pragma once

#include <engine/math/types.h>
#include <engine/util/functions.h>

// Misc
static inline f32_t min(f32_t a, f32_t b) {
	return a < b ? a : b;
}

static inline f32_t max(f32_t a, f32_t b) {
	return a > b ? a : b;
}

static inline f32_t to_rad(f32_t deg) {
	return deg * PI / 180.0f;
}

static inline f32_t to_deg(f32_t rad) {
	return rad * 180.0f / PI;
}

// vec2_t
static inline vec2_t vec2_init_f(f32_t f) {
	const vec2_t v = {f, f};
	return v;
}

static inline vec2_t vec2_init_2f(f32_t x, f32_t y) {
	const vec2_t v = {x, y};
	return v;
}

static inline vec2_t vec2_mix(vec2_t a, vec2_t b, vec2_t f) {
	return a * (vec2_init_f(1.0f) - f) + b * f;
}

// vec3_t
static inline vec3_t vec3_init_f(f32_t f) {
	vec3_t v = {f, f, f};
	return v;
}

static inline vec3_t vec3_init_3f(f32_t x, f32_t y, f32_t z) {
	vec3_t v = {x, y, z};
	return v;
}

static inline f32_t vec3_dot(vec3_t a, vec3_t b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline f32_t vec3_length(vec3_t a) {
	const f32_t len_sq = vec3_dot(a, a);
	return sqrtf(len_sq);
}

static inline vec3_t vec3_normalize(vec3_t a) {
	const f32_t  len = vec3_length(a);
	const vec3_t y   = vec3_init_f(len);
	return a / y;
}

static inline vec3_t vec3_cross(vec3_t a, vec3_t b) {
	vec3_t dest;
	dest.x = a.y * b.z - a.z * b.y;
	dest.y = a.z * b.x - a.x * b.z;
	dest.z = a.x * b.y - a.y * b.x;
	return dest;
}

static inline vec3_t vec3_mix(vec3_t a, vec3_t b, f32_t f) {
	return a * vec3_init_f(1.0f - f) + b * vec3_init_f(f);
}

static inline vec3_t vec3_clamp(vec3_t a, vec3_t min_v, vec3_t max_v) {
	vec3_t out;
	out.x = max(min_v.x, min(a.x, max_v.x));
	out.y = max(min_v.y, min(a.y, max_v.y));
	out.z = max(min_v.z, min(a.z, max_v.z));
	return out;
}

// vec4_t
static inline vec4_t vec4_init_f(f32_t f) {
	const vec4_t v = {f, f, f, f};
	return v;
}

static inline vec4_t vec4_init_4f(f32_t x, f32_t y, f32_t z, f32_t w) {
	const vec4_t v = {x, y, z, w};
	return v;
}

static inline f32_t vec4_dot(vec4_t a, vec4_t b) {
	return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

// mat3_t
mat3_t mat3_init_mat4(const mat4_t *m4);
vec3_t mat3_mul_vec3(const mat3_t *m, vec3_t v);

// mat4_t
mat4_t mat4_init_f(f32_t f);
mat4_t mat4_mul_mat4(const mat4_t *a, const mat4_t *b);
vec4_t mat4_mul_vec4(const mat4_t *a, vec4_t b);
mat4_t mat4_invert(const mat4_t *mat);
mat4_t look_at(vec3_t eye, vec3_t center, vec3_t up);
mat4_t perspective(f32_t fovy, f32_t aspect, f32_t near, f32_t far);

// Sphere
static inline sphere_t sphere_init(vec3_t center, f32_t radius) {
	const sphere_t s = {center, radius};
	return s;
}

// Plane
static inline plane_t plane_init(vec3_t normal, f32_t offset) {
	const plane_t p = {normal, offset};
	return p;
}

static inline f32_t plane_point_distance(const plane_t *plane, vec3_t point) {
	return vec3_dot(plane->normal, point) - plane->offset;
}

// Ray intersection
bool_t ray_cast_sphere(const ray_t *ray, const sphere_t *sphere, f32_t t_min, f32_t t_max, ray_hit_t *hit);
bool_t ray_cast_plane(const ray_t *ray, const plane_t *plane, f32_t t_min, f32_t t_max, ray_hit_t *hit);

// Other
static inline vec3_t random_in_init_sphere() {
	vec3_t out;
	do {
		out = vec3_init_3f(rand_0f_to_1f(), rand_0f_to_1f(), rand_0f_to_1f());
		out = out * 2.0f - 1.0f;
	} while(vec3_length(out) <= 1.0f);
	return out;
}
