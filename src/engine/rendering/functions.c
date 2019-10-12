#include <engine/rendering/functions.h>

static f32_t schlick(f32_t cosine, f32_t ref_idx) {
	f32_t r0 = (1.0f - ref_idx) / (1.0f + ref_idx);
	r0 *= r0;
	return r0 + (1.0f - r0) * powf(1.0f - cosine, 5.0f);
}

bool_t lambertian_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                          ray_t *scattered_ray, f32_t *pdf) {
	(void)in_ray;
	const vec3_t target = hit->normal + random_in_unit_sphere();
	*scattered_ray = ray_init(hit->point, vec3_normalize(target));
	*attenuation = mtl->albedo_texture.callback(&mtl->albedo_texture, hit->uv, hit->point);
	*pdf = vec3_dot(hit->normal, scattered_ray->direction) / PI;
	return TRUE;
}

bool_t metal_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                     ray_t *scattered_ray, f32_t *pdf) {
	const vec3_t reflected = vec3_normalize(vec3_reflect(in_ray->direction, hit->normal));
	const f32_t metal_fuzz =
	    mtl->metal_fuzz_texture.callback(&mtl->metal_fuzz_texture, vec2_init_f(0.0f), hit->point).x;
	const vec3_t target = reflected + random_in_unit_sphere() * metal_fuzz;
	*scattered_ray = ray_init(hit->point, vec3_normalize(target));
	*attenuation = mtl->albedo_texture.callback(&mtl->albedo_texture, vec2_init_f(0.0f), hit->point);
	*pdf = 0.0f; // TODO
	return vec3_dot(reflected, hit->normal) > 0.0f;
}

bool_t dielectric_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                          ray_t *scattered_ray, f32_t *pdf) {
	vec3_t outward_normal;
	const vec3_t reflected = vec3_normalize(vec3_reflect(in_ray->direction, hit->normal));
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

	*pdf = 0.0f; // TODO
	return TRUE;
}

aabb_t render_queue_compute_aabb(const render_queue_t *rqueue) {
	aabb_t out = aabb_init(vec3_init_f(0.0f), vec3_init_f(0.0001f));

	if(rqueue->renderable_count > 0) {
		out = sphere_compute_aabb(&rqueue->renderables[0].shape.sphere);

		for(u32_t i = 1; i < rqueue->renderable_count; ++i) {
			const aabb_t box = sphere_compute_aabb(&rqueue->renderables[i].shape.sphere);
			out = aabb_union(&out, &box);
		}
	}

	return out;
}

bool_t render_queue_closest_hit(const render_queue_t *rgraph, const ray_t *ray, ray_hit_t *closest_hit,
                                const material_t **hit_mtl) {
	closest_hit->t = INFINITY;
	bool_t has_hit = FALSE;
	const f32_t t_min = 0.001;
	const f32_t t_max = 5000.0f;

	for(u32_t i = 0; i < rgraph->renderable_count; ++i) {
		const renderable_t *renderable = &rgraph->renderables[i];
		ray_hit_t hit;

		if(renderable->shape_type == RENDERABLE_SHAPE_TYPE_SPHERE) {
			const sphere_t *sphere = &renderable->shape.sphere;

			sphere_t new_sphere = sphere_init(sphere->center + renderable->world_transform.translation, sphere->radius);
			const sphere_t old_sphere =
			    sphere_init(sphere->center + renderable->previous_world_transform.translation, sphere->radius);

			if(!vec3_eq(new_sphere.center, old_sphere.center)) {
				new_sphere.center = vec3_mix(new_sphere.center, old_sphere.center, rand_0f_to_1f());
			}
			if(ray_cast_sphere(ray, &new_sphere, t_min, t_max, &hit) && hit.t < closest_hit->t) {
				*closest_hit = hit;
				*hit_mtl = &renderable->material;
				has_hit = TRUE;
			}
		} else if(renderable->shape_type == RENDERABLE_SHAPE_TYPE_MESH) {
			// TODO won't work with scale
			const mesh_t *mesh = &renderable->shape.mesh;
			bool_t mesh_has_hit = FALSE;

			const mat4_t trf = mat4_init_transform(&renderable->world_transform);
			const mat4_t inv_trf = mat4_invert(&trf);
			const mat3_t rot = mat3_init_mat4(&trf);
			const mat3_t inv_rot = mat3_init_mat4(&inv_trf);
			ray_t new_ray;
			new_ray.origin =
			    mat4_mul_vec4(&inv_trf, vec4_init_4f(ray->origin.x, ray->origin.y, ray->origin.z, 1.0f)).xyz;
			new_ray.direction = mat3_mul_vec3(&inv_rot, ray->direction);

			for(u32_t t = 0; t < mesh->triangle_count; ++t) {
				if(ray_cast_triangle(&new_ray, &mesh->triangles[t], t_min, t_max, &hit) && hit.t < closest_hit->t) {
					*closest_hit = hit;
					*hit_mtl = &renderable->material;
					mesh_has_hit = TRUE;
				}
			}

			for(u32_t t = 0; t < mesh->quad_count; ++t) {
				if(ray_cast_quad(&new_ray, &mesh->quads[t], t_min, t_max, &hit) && hit.t < closest_hit->t) {
					*closest_hit = hit;
					*hit_mtl = &renderable->material;
					mesh_has_hit = TRUE;
				}
			}

			if(mesh_has_hit) {
				closest_hit->normal = mat3_mul_vec3(&rot, closest_hit->normal);
				closest_hit->point = mat4_mul_vec4(&trf, vec4_init_vec3(closest_hit->point, 1.0f)).xyz;
				has_hit = TRUE;
			}
		} else {
			assert(0);
		}
	}

	return has_hit;
}
