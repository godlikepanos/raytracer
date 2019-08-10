#include <engine/renderer/rendering_functions.h>

static f32_t schlick(f32_t cosine, f32_t ref_idx) {
	f32_t r0 = (1.0f - ref_idx) / (1.0f + ref_idx);
	r0 *= r0;
	return r0 + (1.0f - r0) * powf(1.0f - cosine, 5.0f);
}

bool_t lambertian_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                          ray_t *scattered_ray) {
	(void)in_ray;
	const vec3_t target = hit->normal + random_in_unit_sphere();
	*scattered_ray = ray_init(hit->point, vec3_normalize(target));
	*attenuation = mtl->albedo_texture.callback(&mtl->albedo_texture, vec2_init_f(0.0f), hit->point);
	return TRUE;
}

bool_t metal_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                     ray_t *scattered_ray) {
	const vec3_t reflected = vec3_normalize(vec3_reflect(in_ray->direction, hit->normal));
	const f32_t metal_fuzz =
	    mtl->metal_fuzz_texture.callback(&mtl->metal_fuzz_texture, vec2_init_f(0.0f), hit->point).x;
	const vec3_t target = reflected + random_in_unit_sphere() * metal_fuzz;
	*scattered_ray = ray_init(hit->point, vec3_normalize(target));
	*attenuation = mtl->albedo_texture.callback(&mtl->albedo_texture, vec2_init_f(0.0f), hit->point);
	return vec3_dot(reflected, hit->normal) > 0.0f;
}

bool_t dielectric_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                          ray_t *scattered_ray) {
	vec3_t outward_normal;
	const vec3_t reflected = vec3_reflect(in_ray->direction, hit->normal);
	f32_t ni_over_nt;
	f32_t cosine;
	*attenuation = vec3_init_3f(1.0f, 1.0f, 1.0f);
	const f32_t refl_idx =
	    mtl->dielectric_reflection_index.callback(&mtl->dielectric_reflection_index, vec2_init_f(0.0f), hit->point).x;
	if(vec3_dot(in_ray->direction, hit->normal) > 0.0f) {
		outward_normal = -hit->normal;
		ni_over_nt = refl_idx;
		cosine = refl_idx * vec3_dot(in_ray->direction, hit->normal);
	} else {
		outward_normal = hit->normal;
		ni_over_nt = 1.0f / refl_idx;
		cosine = -vec3_dot(in_ray->direction, hit->normal);
	}

	vec3_t refracted;
	f32_t reflect_probe;
	if(vec3_refract(in_ray->direction, outward_normal, ni_over_nt, &refracted)) {
		reflect_probe = schlick(cosine, refl_idx);
	} else {
		reflect_probe = 1.0f;
	}

	if(rand_0f_to_1f() < reflect_probe) {
		*scattered_ray = ray_init(hit->point, reflected);
	} else {
		*scattered_ray = ray_init(hit->point, refracted);
	}

	return TRUE;
}

aabb_t render_queue_compute_aabb(const render_queue_t *rqueue) {
	aabb_t out = aabb_init(vec3_init_f(0.0f), vec3_init_f(0.0001f));

	if(rqueue->sphere_count > 0) {
		out = sphere_compute_aabb(&rqueue->spheres[0].shape.sphere);

		for(u32_t i = 1; i < rqueue->sphere_count; ++i) {
			const aabb_t box = sphere_compute_aabb(&rqueue->spheres[i].shape.sphere);
			out = aabb_union(&out, &box);
		}
	}

	return out;
}
