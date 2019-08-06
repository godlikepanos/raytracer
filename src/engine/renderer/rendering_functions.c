#include <engine/renderer/rendering_functions.h>

bool_t lambertian_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                          ray_t *scattered_ray) {
	(void)in_ray;
	const vec3_t target = hit->normal + random_in_unit_sphere();
	*scattered_ray = ray_init(hit->point, vec3_normalize(target));
	*attenuation = mtl->albedo;
	return TRUE;
}

bool_t metal_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                     ray_t *scattered_ray) {
	const vec3_t reflected = vec3_normalize(vec3_reflect(in_ray->direction, hit->normal));
	const vec3_t target = reflected + random_in_unit_sphere() * mtl->metal_fuzz;
	*scattered_ray = ray_init(hit->point, vec3_normalize(target));
	*attenuation = mtl->albedo;
	return vec3_dot(reflected, hit->normal) > 0.0f;
}

bool_t dielectric_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                          ray_t *scattered_ray) {
	vec3_t outward_normal;
	const vec3_t reflected = vec3_reflect(in_ray->direction, hit->normal);
	f32_t ni_over_nt;
	f32_t cosine;
	*attenuation = vec3_init_3f(1.0f, 1.0f, 1.0f);
	if(vec3_dot(in_ray->direction, hit->normal) > 0.0f) {
		outward_normal = -hit->normal;
		ni_over_nt = mtl->dielectric_refl_idx;
		cosine = mtl->dielectric_refl_idx * vec3_dot(in_ray->direction, hit->normal);
	} else {
		outward_normal = hit->normal;
		ni_over_nt = 1.0f / mtl->dielectric_refl_idx;
		cosine = -vec3_dot(in_ray->direction, hit->normal);
	}

	vec3_t refracted;
	f32_t reflect_probe;
	if(vec3_refract(in_ray->direction, outward_normal, ni_over_nt, &refracted)) {
		reflect_probe = schlick(cosine, mtl->dielectric_refl_idx);
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
