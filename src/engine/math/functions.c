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
	f32_t tmp[12];
	mat4_t m4;

	tmp[0] = mat->m[2][2] * mat->m[3][3];
	tmp[1] = mat->m[3][2] * mat->m[2][3];
	tmp[2] = mat->m[1][2] * mat->m[3][3];
	tmp[3] = mat->m[3][2] * mat->m[1][3];
	tmp[4] = mat->m[1][2] * mat->m[2][3];
	tmp[5] = mat->m[2][2] * mat->m[1][3];
	tmp[6] = mat->m[0][2] * mat->m[3][3];
	tmp[7] = mat->m[3][2] * mat->m[0][3];
	tmp[8] = mat->m[0][2] * mat->m[2][3];
	tmp[9] = mat->m[2][2] * mat->m[0][3];
	tmp[10] = mat->m[0][2] * mat->m[1][3];
	tmp[11] = mat->m[1][2] * mat->m[0][3];

	m4.m[0][0] = tmp[0] * mat->m[1][1] + tmp[3] * mat->m[2][1] + tmp[4] * mat->m[3][1];
	m4.m[0][0] -= tmp[1] * mat->m[1][1] + tmp[2] * mat->m[2][1] + tmp[5] * mat->m[3][1];
	m4.m[0][1] = tmp[1] * mat->m[0][1] + tmp[6] * mat->m[2][1] + tmp[9] * mat->m[3][1];
	m4.m[0][1] -= tmp[0] * mat->m[0][1] + tmp[7] * mat->m[2][1] + tmp[8] * mat->m[3][1];
	m4.m[0][2] = tmp[2] * mat->m[0][1] + tmp[7] * mat->m[1][1] + tmp[10] * mat->m[3][1];
	m4.m[0][2] -= tmp[3] * mat->m[0][1] + tmp[6] * mat->m[1][1] + tmp[11] * mat->m[3][1];
	m4.m[0][3] = tmp[5] * mat->m[0][1] + tmp[8] * mat->m[1][1] + tmp[11] * mat->m[2][1];
	m4.m[0][3] -= tmp[4] * mat->m[0][1] + tmp[9] * mat->m[1][1] + tmp[10] * mat->m[2][1];
	m4.m[1][0] = tmp[1] * mat->m[1][0] + tmp[2] * mat->m[2][0] + tmp[5] * mat->m[3][0];
	m4.m[1][0] -= tmp[0] * mat->m[1][0] + tmp[3] * mat->m[2][0] + tmp[4] * mat->m[3][0];
	m4.m[1][1] = tmp[0] * mat->m[0][0] + tmp[7] * mat->m[2][0] + tmp[8] * mat->m[3][0];
	m4.m[1][1] -= tmp[1] * mat->m[0][0] + tmp[6] * mat->m[2][0] + tmp[9] * mat->m[3][0];
	m4.m[1][2] = tmp[3] * mat->m[0][0] + tmp[6] * mat->m[1][0] + tmp[11] * mat->m[3][0];
	m4.m[1][2] -= tmp[2] * mat->m[0][0] + tmp[7] * mat->m[1][0] + tmp[10] * mat->m[3][0];
	m4.m[1][3] = tmp[4] * mat->m[0][0] + tmp[9] * mat->m[1][0] + tmp[10] * mat->m[2][0];
	m4.m[1][3] -= tmp[5] * mat->m[0][0] + tmp[8] * mat->m[1][0] + tmp[11] * mat->m[2][0];

	tmp[0] = mat->m[2][0] * mat->m[3][1];
	tmp[1] = mat->m[3][0] * mat->m[2][1];
	tmp[2] = mat->m[1][0] * mat->m[3][1];
	tmp[3] = mat->m[3][0] * mat->m[1][1];
	tmp[4] = mat->m[1][0] * mat->m[2][1];
	tmp[5] = mat->m[2][0] * mat->m[1][1];
	tmp[6] = mat->m[0][0] * mat->m[3][1];
	tmp[7] = mat->m[3][0] * mat->m[0][1];
	tmp[8] = mat->m[0][0] * mat->m[2][1];
	tmp[9] = mat->m[2][0] * mat->m[0][1];
	tmp[10] = mat->m[0][0] * mat->m[1][1];
	tmp[11] = mat->m[1][0] * mat->m[0][1];

	m4.m[2][0] = tmp[0] * mat->m[1][3] + tmp[3] * mat->m[2][3] + tmp[4] * mat->m[3][3];
	m4.m[2][0] -= tmp[1] * mat->m[1][3] + tmp[2] * mat->m[2][3] + tmp[5] * mat->m[3][3];
	m4.m[2][1] = tmp[1] * mat->m[0][3] + tmp[6] * mat->m[2][3] + tmp[9] * mat->m[3][3];
	m4.m[2][1] -= tmp[0] * mat->m[0][3] + tmp[7] * mat->m[2][3] + tmp[8] * mat->m[3][3];
	m4.m[2][2] = tmp[2] * mat->m[0][3] + tmp[7] * mat->m[1][3] + tmp[10] * mat->m[3][3];
	m4.m[2][2] -= tmp[3] * mat->m[0][3] + tmp[6] * mat->m[1][3] + tmp[11] * mat->m[3][3];
	m4.m[2][3] = tmp[5] * mat->m[0][3] + tmp[8] * mat->m[1][3] + tmp[11] * mat->m[2][3];
	m4.m[2][3] -= tmp[4] * mat->m[0][3] + tmp[9] * mat->m[1][3] + tmp[10] * mat->m[2][3];
	m4.m[3][0] = tmp[2] * mat->m[2][2] + tmp[5] * mat->m[3][2] + tmp[1] * mat->m[1][2];
	m4.m[3][0] -= tmp[4] * mat->m[3][2] + tmp[0] * mat->m[1][2] + tmp[3] * mat->m[2][2];
	m4.m[3][1] = tmp[8] * mat->m[3][2] + tmp[0] * mat->m[0][2] + tmp[7] * mat->m[2][2];
	m4.m[3][1] -= tmp[6] * mat->m[2][2] + tmp[9] * mat->m[3][2] + tmp[1] * mat->m[0][2];
	m4.m[3][2] = tmp[6] * mat->m[1][2] + tmp[11] * mat->m[3][2] + tmp[3] * mat->m[0][2];
	m4.m[3][2] -= tmp[10] * mat->m[3][2] + tmp[2] * mat->m[0][2] + tmp[7] * mat->m[1][2];
	m4.m[3][3] = tmp[10] * mat->m[2][2] + tmp[4] * mat->m[0][2] + tmp[9] * mat->m[1][2];
	m4.m[3][3] -= tmp[8] * mat->m[1][2] + tmp[11] * mat->m[2][2] + tmp[5] * mat->m[0][2];

	f32_t det =
	    mat->m[0][0] * m4.m[0][0] + mat->m[1][0] * m4.m[0][1] + mat->m[2][0] * m4.m[0][2] + mat->m[3][0] * m4.m[0][3];
	det = 1.0 / det;
	m4.m[0] *= det;
	m4.m[1] *= det;
	m4.m[2] *= det;
	m4.m[3] *= det;
	return m4;
}

mat4_t look_at(vec3_t eye, vec3_t center, vec3_t up) {
	const vec3_t vdir = vec3_normalize(center - eye);
	const vec3_t vside = vec3_normalize(vec3_cross(vdir, up));
	const vec3_t vup = vec3_normalize(vec3_cross(vside, vdir));
	const vec3_t x = vside;
	const vec3_t y = vup;
	const vec3_t z = -vdir;

	mat4_t dest;
	dest.m[0][0] = x.x;
	dest.m[0][1] = y.x;
	dest.m[0][2] = z.x;
	dest.m[0][3] = eye.x;
	dest.m[1][0] = x.y;
	dest.m[1][1] = y.y;
	dest.m[1][2] = z.y;
	dest.m[1][3] = eye.y;
	dest.m[2][0] = x.z;
	dest.m[2][1] = y.z;
	dest.m[2][2] = z.z;
	dest.m[2][3] = eye.z;
	dest.m[3][0] = 0.0f;
	dest.m[3][1] = 0.0f;
	dest.m[3][2] = 0.0f;
	dest.m[3][3] = 1.0f;

	return mat4_invert(&dest);
}

mat4_t perspective(f32_t fovy, f32_t aspect, f32_t near, f32_t far) {
	const f32_t g = near - far;
	const f32_t f = 1.0f / tanf(fovy / 2.0f);

	mat4_t dest;
	dest.m[0][0] = f / aspect;
	dest.m[0][1] = 0.0f;
	dest.m[0][2] = 0.0f;
	dest.m[0][3] = 0.0f;
	dest.m[1][0] = 0.0f;
	dest.m[1][1] = f;
	dest.m[1][2] = 0.0f;
	dest.m[1][3] = 0.0f;
	dest.m[2][0] = 0.0f;
	dest.m[2][1] = 0.0f;
	dest.m[2][2] = far / g;
	dest.m[2][3] = (far * near) / g;
	dest.m[3][0] = 0.0f;
	dest.m[3][1] = 0.0f;
	dest.m[3][2] = -1.0f;
	dest.m[3][3] = 0.0f;

	return dest;
}

bool_t ray_cast_sphere(const ray_t *ray, const sphere_t *sphere, f32_t t_min, f32_t t_max, ray_hit_t *hit) {
	vec3_t oc = ray->origin - sphere->center;
	f32_t a = vec3_dot(ray->direction, ray->direction);
	f32_t b = vec3_dot(oc, ray->direction);
	f32_t c = vec3_dot(oc, oc) - sphere->radius * sphere->radius;
	f32_t d = b * b - a * c;

	if(d > 0.0f) {
		f32_t tmp = (-b - sqrtf(d)) / a;
		if(tmp < t_max && tmp > t_min) {
			hit->t = tmp;
			hit->point = ray->origin + ray->direction * hit->t;
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
			hit->point = ray->origin + ray->direction * hit->t;
			hit->normal = plane->normal;
			return TRUE;
		}
	}

	return FALSE;
}

bool_t ray_intersects_aabb(const ray_t *ray, const aabb_t *box, f32_t tmin, f32_t tmax) {
	for(u32_t a = 0; a < 3; ++a) {
		const f32_t l = (box->min[a] - ray->origin[a]) / ray->direction[a];
		const f32_t r = (box->max[a] - ray->origin[a]) / ray->direction[a];
		const f32_t t0 = MIN(l, r);
		const f32_t t1 = MAX(l, r);
		tmin = MAX(t0, tmin);
		tmax = MIN(t1, tmax);
		if(tmax <= tmin) {
			return FALSE;
		}
	}
	return TRUE;
}
