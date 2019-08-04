#include <engine/math/functions.h>

mat3_t mat3_init_mat4(const mat4_t *m4) {
	mat3_t out;
	for(u32_t j = 0; j < 3; ++j) {
		for(u32_t i = 0; i < 3; ++i) {
			out.m[j][i] = m4->m[j][i];
		}
	}
	return out;
}

vec3_t mat3_mul_vec3(const mat3_t *m, vec3_t v) {
	vec3_t out;
	for(u32_t j = 0; j < 3; j++) {
		f32_t sum = 0.0f;
		for(u32_t i = 0; i < 3; i++) {
			sum += m->m[j][i] * v[i];
		}
		out[j] = sum;
	}
	return out;
}

mat4_t mat4_init_f(f32_t f) {
	mat4_t m;
	for(u32_t i = 0; i < 4; ++i)
		m.m[i] = vec4_init_f(f);
	return m;
}

mat4_t mat4_mul_mat4(const mat4_t *a, const mat4_t *b) {
	mat4_t m;
	for(u32_t j = 0; j < 4; j++) {
		for(u32_t i = 0; i < 4; i++) {
			m.m[j][i] = 0.0f;
			for(u32_t k = 0; k < 4; k++) {
				m.m[j][i] += a->m[j][k] * b->m[k][i];
			}
		}
	}

	return m;
}

vec4_t mat4_mul_vec4(const mat4_t *m, vec4_t v) {
	vec4_t out;
	for(u32_t j = 0; j < 4; j++) {
		f32_t sum = 0.0f;
		for(u32_t i = 0; i < 4; i++) {
			sum += m->m[j][i] * v[i];
		}
		out[j] = sum;
	}
	return out;
}

mat4_t mat4_invert(const mat4_t *mat) {
	f32_t t[6];
	f32_t a = mat->m[0][0], b = mat->m[0][1], c = mat->m[0][2], d = mat->m[0][3], e = mat->m[1][0], f = mat->m[1][1],
	      g = mat->m[1][2], h = mat->m[1][3], i = mat->m[2][0], j = mat->m[2][1], k = mat->m[2][2], l = mat->m[2][3],
	      m = mat->m[3][0], n = mat->m[3][1], o = mat->m[3][2], p = mat->m[3][3];

	t[0] = k * p - o * l;
	t[1] = j * p - n * l;
	t[2] = j * o - n * k;
	t[3] = i * p - m * l;
	t[4] = i * o - m * k;
	t[5] = i * n - m * j;

	mat4_t dest;
	dest.m[0][0] = f * t[0] - g * t[1] + h * t[2];
	dest.m[1][0] = -(e * t[0] - g * t[3] + h * t[4]);
	dest.m[2][0] = e * t[1] - f * t[3] + h * t[5];
	dest.m[3][0] = -(e * t[2] - f * t[4] + g * t[5]);

	dest.m[0][1] = -(b * t[0] - c * t[1] + d * t[2]);
	dest.m[1][1] = a * t[0] - c * t[3] + d * t[4];
	dest.m[2][1] = -(a * t[1] - b * t[3] + d * t[5]);
	dest.m[3][1] = a * t[2] - b * t[4] + c * t[5];

	t[0] = g * p - o * h;
	t[1] = f * p - n * h;
	t[2] = f * o - n * g;
	t[3] = e * p - m * h;
	t[4] = e * o - m * g;
	t[5] = e * n - m * f;

	dest.m[0][2] = b * t[0] - c * t[1] + d * t[2];
	dest.m[1][2] = -(a * t[0] - c * t[3] + d * t[4]);
	dest.m[2][2] = a * t[1] - b * t[3] + d * t[5];
	dest.m[3][2] = -(a * t[2] - b * t[4] + c * t[5]);

	t[0] = g * l - k * h;
	t[1] = f * l - j * h;
	t[2] = f * k - j * g;
	t[3] = e * l - i * h;
	t[4] = e * k - i * g;
	t[5] = e * j - i * f;

	dest.m[0][3] = -(b * t[0] - c * t[1] + d * t[2]);
	dest.m[1][3] = a * t[0] - c * t[3] + d * t[4];
	dest.m[2][3] = -(a * t[1] - b * t[3] + d * t[5]);
	dest.m[3][3] = a * t[2] - b * t[4] + c * t[5];

	const f32_t det = 1.0f / (a * dest.m[0][0] + b * dest.m[1][0] + c * dest.m[2][0] + d * dest.m[3][0]);

	dest.m[0] *= vec4_init_f(det);
	dest.m[1] *= vec4_init_f(det);
	dest.m[2] *= vec4_init_f(det);
	dest.m[3] *= vec4_init_f(det);

	return dest;
}

mat4_t look_at(vec3_t eye, vec3_t center, vec3_t up) {
	const vec3_t f = vec3_normalize(center - eye);
	const vec3_t s = vec3_normalize(vec3_cross(f, up));
	const vec3_t u = vec3_cross(s, f);

	mat4_t dest;
	dest.m[0][0] = s.x;
	dest.m[0][1] = u.x;
	dest.m[0][2] = -f.x;
	dest.m[1][0] = s.y;
	dest.m[1][1] = u.y;
	dest.m[1][2] = -f.y;
	dest.m[2][0] = s.z;
	dest.m[2][1] = u.z;
	dest.m[2][2] = -f.z;
	dest.m[3][0] = -vec3_dot(s, eye);
	dest.m[3][1] = -vec3_dot(u, eye);
	dest.m[3][2] = vec3_dot(f, eye);
	dest.m[0][3] = 0.0f;
	dest.m[1][3] = 0.0f;
	dest.m[2][3] = 0.0f;
	dest.m[3][3] = 1.0f;

	return dest;
}

mat4_t perspective(f32_t fovy, f32_t aspect, f32_t near, f32_t far) {
	const f32_t f  = 1.0f / tanf(fovy * 0.5f);
	const f32_t fn = 1.0f / (near - far);

	mat4_t dest  = mat4_init_f(0.0f);
	dest.m[0][0] = f / aspect;
	dest.m[1][1] = f;
	dest.m[2][2] = (near + far) * fn;
	dest.m[2][3] = -1.0f;
	dest.m[3][2] = 2.0f * near * far * fn;

	return dest;
}

bool_t ray_cast_sphere(const ray_t *ray, const sphere_t *sphere, f32_t t_min, f32_t t_max, ray_hit_t *hit) {
	vec3_t oc = ray->origin - sphere->center;
	f32_t  a  = vec3_dot(ray->direction, ray->direction);
	f32_t  b  = vec3_dot(oc, ray->direction);
	f32_t  c  = vec3_dot(oc, oc) - sphere->radius * sphere->radius;
	f32_t  d  = b * b - a * c;

	if(d > 0.0f) {
		f32_t tmp = (-b - sqrtf(d)) / a;
		if(tmp < t_max && tmp > t_min) {
			hit->t      = tmp;
			hit->point  = ray->origin + ray->direction * hit->t;
			hit->normal = (hit->point - sphere->center) / sphere->radius;
			return TRUE;
		}
	}

	return FALSE;
}

bool_t ray_cast_plane(const ray_t *ray, const plane_t *plane, f32_t t_min, f32_t t_max, ray_hit_t *hit) {
	const f32_t d = plane_point_distance(plane, ray->origin);
	const f32_t a = vec3_dot(plane->normal, ray->direction);

	if(d > 0.0f && a < 0.0f) {
		const f32_t tmp = -d / a;
		if(tmp < t_max && tmp > t_min) {
			hit->t = -d / a;
			assert(hit->t > 0.0f);
			hit->point  = ray->origin + ray->direction * hit->t;
			hit->normal = plane->normal;
			return TRUE;
		}
	}

	return FALSE;
}
