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

static inline bool_t vec3_eq(vec3_t a, vec3_t b) {
	return a.x == b.x && a.y == b.y && a.z == b.z;
}

static inline vec3_t vec3_min(vec3_t a, vec3_t b) {
	vec3_t c;
	c.x = MIN(a.x, b.x);
	c.y = MIN(a.y, b.y);
	c.z = MIN(a.z, b.z);
	return c;
}

static inline vec3_t vec3_max(vec3_t a, vec3_t b) {
	vec3_t c;
	c.x = MAX(a.x, b.x);
	c.y = MAX(a.y, b.y);
	c.z = MAX(a.z, b.z);
	return c;
}

static inline f32_t vec3_dot(vec3_t a, vec3_t b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline f32_t vec3_length_squared(vec3_t a) {
	const f32_t len_sq = vec3_dot(a, a);
	return len_sq;
}

static inline f32_t vec3_length(vec3_t a) {
	return sqrtf(vec3_length_squared(a));
}

static inline vec3_t vec3_normalize(vec3_t a) {
	const f32_t len = vec3_length(a);
	return a / len;
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

static inline vec3_t vec3_reflect(vec3_t v, vec3_t n) {
	return v - 2.0f * vec3_dot(v, n) * n;
}

static inline bool_t vec3_refract(vec3_t v, vec3_t n, f32_t ni_over_nt, vec3_t *refracted) {
	const vec3_t unit_v = vec3_normalize(v);
	const f32_t dt = vec3_dot(unit_v, n);
	const f32_t discriminant = 1.0f - ni_over_nt * ni_over_nt * (1.0f - dt * dt);
	if(discriminant > 0.0f) {
		*refracted = ni_over_nt * (unit_v - n * dt) - n * sqrtf(discriminant);
		return TRUE;
	} else {
		return FALSE;
	}
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

static inline vec4_t vec4_init_vec3(vec3_t v3, f32_t w) {
	const vec4_t v = {v3.x, v3.y, v3.z, w};
	return v;
}

static inline f32_t vec4_dot(vec4_t a, vec4_t b) {
	return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

// mat3_t
mat3_t mat3_init_identity();
mat3_t mat3_init_mat4(const mat4_t *m4);
mat3_t mat3_init_axis_angles(vec3_t axis, f32_t angle);
mat3_t mat3_init_columns(vec3_t x, vec3_t y, vec3_t z);
vec3_t mat3_mul_vec3(const mat3_t *m, vec3_t v);
mat3_t mat3_orthonormal_basis(vec3_t z);
vec3_t mat3_get_column(const mat3_t *m, u32_t column);

// mat4_t
mat4_t mat4_init_f(f32_t f);
mat4_t mat4_init_transform(const transform_t *trf);
mat4_t mat4_mul_mat4(const mat4_t *a, const mat4_t *b);
vec4_t mat4_mul_vec4(const mat4_t *a, vec4_t b);
mat4_t mat4_invert(const mat4_t *mat);
mat4_t look_at(vec3_t eye, vec3_t center, vec3_t up);
mat4_t perspective(f32_t fovy, f32_t aspect, f32_t near, f32_t far);

// transform_t
static inline transform_t transform_init_trs(vec3_t translation, const mat3_t *rotation, vec3_t scale) {
	const transform_t out = {translation, *rotation, scale};
	return out;
}

static inline transform_t transform_init_t(vec3_t translation) {
	const mat3_t identity = mat3_init_identity();
	return transform_init_trs(translation, &identity, vec3_init_f(1.0f));
}

static inline transform_t transform_init_identity() {
	return transform_init_t(vec3_init_f(0.0f));
}

// sphere_t
static inline sphere_t sphere_init(vec3_t center, f32_t radius) {
	const sphere_t s = {center, radius};
	return s;
}

static inline aabb_t sphere_compute_aabb(const sphere_t *s) {
	aabb_t box;
	box.min = s->center - s->radius;
	box.max = s->center + s->radius;
	return box;
}

static inline vec2_t sphere_compute_uv(vec3_t point_in_sphere) {
	const f32_t phi = atan2(point_in_sphere.z, point_in_sphere.x);
	const f32_t theta = asin(point_in_sphere.y);
	vec2_t uv;
	uv.x = 1.0f - (phi + PI) / (2.0f * PI);
	uv.y = (theta + PI / 2.0f) / PI;
	return uv;
}

// plane_t
static inline plane_t plane_init(vec3_t normal, f32_t offset) {
	const plane_t p = {normal, offset};
	return p;
}

static inline f32_t plane_point_distance(const plane_t *plane, vec3_t point) {
	return vec3_dot(plane->normal, point) - plane->offset;
}

// ray_t
static inline ray_t ray_init(vec3_t point, vec3_t normal) {
	const f32_t normal_len = vec3_length(normal);
	assert(normal_len <= 1.001f);
	(void)normal_len;
	const ray_t out = {point, normal};
	return out;
}

bool_t ray_cast_sphere(const ray_t *ray, const sphere_t *sphere, f32_t t_min, f32_t t_max, ray_hit_t *hit);
bool_t ray_cast_plane(const ray_t *ray, const plane_t *plane, f32_t t_min, f32_t t_max, ray_hit_t *hit);
bool_t ray_cast_triangle(const ray_t *ray, const triangle_t *tri, f32_t t_min, f32_t t_max, ray_hit_t *hit);
bool_t ray_cast_quad(const ray_t *ray, const quadrilateral_t *quad, f32_t t_min, f32_t t_max, ray_hit_t *hit);
bool_t ray_cast_aabb(const ray_t *ray, const aabb_t *box, f32_t t_min, f32_t t_max, ray_hit_t *hit);

bool_t ray_intersects_aabb(const ray_t *ray, const aabb_t *box, f32_t tmin, f32_t tmax);

// aabb_t
static inline aabb_t aabb_init(vec3_t min, vec3_t max) {
	const aabb_t box = {min, max};
	return box;
}

static inline aabb_t aabb_union(const aabb_t *a, const aabb_t *b) {
	aabb_t c;
	c.min = vec3_min(a->min, b->min);
	c.max = vec3_min(a->max, b->max);
	return c;
}

// quadrilateral_t
static inline quadrilateral_t quad_flip(const quadrilateral_t *q) {
	quadrilateral_t out;
	out.vertices[0] = q->vertices[0];
	out.vertices[1] = q->vertices[3];
	out.vertices[2] = q->vertices[2];
	out.vertices[3] = q->vertices[1];
	return out;
}

static inline quadrilateral_t quad_init_xy(vec2_t x, vec2_t y, f32_t z, bool_t flip) {
	quadrilateral_t q;
	q.vertices[0] = vec3_init_3f(x.x, y.x, z);
	q.vertices[1] = vec3_init_3f(x.y, y.x, z);
	q.vertices[2] = vec3_init_3f(x.y, y.y, z);
	q.vertices[3] = vec3_init_3f(x.x, y.y, z);
	return (!flip) ? q : quad_flip(&q);
}

static inline quadrilateral_t quad_init_xz(vec2_t x, f32_t y, vec2_t z, bool_t flip) {
	quadrilateral_t q;
	q.vertices[0] = vec3_init_3f(x.x, y, z.y);
	q.vertices[1] = vec3_init_3f(x.y, y, z.y);
	q.vertices[2] = vec3_init_3f(x.y, y, z.x);
	q.vertices[3] = vec3_init_3f(x.x, y, z.x);
	return (!flip) ? q : quad_flip(&q);
}

static inline quadrilateral_t quad_init_yz(f32_t x, vec2_t y, vec2_t z, bool_t flip) {
	quadrilateral_t q;
	q.vertices[0] = vec3_init_3f(x, y.x, z.y);
	q.vertices[1] = vec3_init_3f(x, y.x, z.x);
	q.vertices[2] = vec3_init_3f(x, y.y, z.x);
	q.vertices[3] = vec3_init_3f(x, y.y, z.y);
	return (!flip) ? q : quad_flip(&q);
}

// Other
static inline vec3_t random_in_unit_sphere() {
	vec3_t out;
	do {
		out = vec3_init_3f(rand_0f_to_1f(), rand_0f_to_1f(), rand_0f_to_1f());
		out = out * 2.0f - 1.0f;
	} while(vec3_length_squared(out) <= 1.0f - EPSILON);
	return out;
}

static inline bool_t f32_is_zero(f32_t num) {
	return fabs(num) <= EPSILON;
}

vec3_t random_cosine_direction();
