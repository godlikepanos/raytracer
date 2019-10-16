#include <engine/rendering/functions.h>

static f32_t schlick(f32_t cosine, f32_t ref_idx) {
	f32_t r0 = (1.0f - ref_idx) / (1.0f + ref_idx);
	r0 *= r0;
	return r0 + (1.0f - r0) * powf(1.0f - cosine, 5.0f);
}

f32_t pdf_compute_value(const pdf_t *pdf, vec3_t direction) {
	assert(f32_is_zero(vec3_length(direction) - 1.0f));

	f32_t out;

	switch(pdf->pdf_type) {
	case PDF_TYPE_COSINE: {
		const f32_t cosine = vec3_dot(direction, mat3_get_column(&pdf->cosine.orthonormal_basis, 2));
		out = (cosine > 0.0f) ? cosine / PI : 0.0f;
		break;
	}
	case PDF_TYPE_HITTABLE: {
		ray_t ray = ray_init(pdf->hittable.ray_origin, direction);
		ray_hit_t hit;
		if(ray_cast_aabb(&ray, &pdf->hittable.box, 0.01f, FLT_MAX, &hit)) {
			const f32_t volume = (pdf->hittable.box.max.x - pdf->hittable.box.min.x)
			                     //* (pdf->hittable.box.max.y - pdf->hittable.box.min.y) TODO
			                     * (pdf->hittable.box.max.z - pdf->hittable.box.min.z);

			const f32_t dist_sq = hit.t * hit.t;
			const f32_t cosine = fabs(vec3_dot(direction, hit.normal));
			out = dist_sq / (cosine * volume);
		} else {
			out = 0.0f;
		}

		break;
	}
	case PDF_TYPE_MIXTURE: {
		const f32_t factor = 1.0f / (f32_t)pdf->mixture.pdf_count;
		out = 0.0f;
		for(u32_t i = 0; i < pdf->mixture.pdf_count; ++i) {
			out += factor * pdf_compute_value(pdf->mixture.pdfs[i], direction);
		}
		break;
	}
	default:
		assert(0);
		out = 0.0f;
	}

	return out;
}

vec3_t pdf_generate(const pdf_t *pdf) {
	vec3_t out = 0.0f;

	switch(pdf->pdf_type) {
	case PDF_TYPE_COSINE:
		out = mat3_mul_vec3(&pdf->cosine.orthonormal_basis, random_cosine_direction());
		break;
	case PDF_TYPE_HITTABLE: {
		const vec3_t random = vec3_mix(pdf->hittable.box.min, pdf->hittable.box.max,
		                               vec3_init_3f(rand_0f_to_1f(), rand_0f_to_1f(), rand_0f_to_1f()));
		out = vec3_normalize(random - pdf->hittable.ray_origin);
		break;
	}
	case PDF_TYPE_MIXTURE: {
		const f32_t r = rand_0f_to_1f();
		const f32_t idxf = r * (f32_t)pdf->mixture.pdf_count;
		const u32_t idx = MIN((u32_t)idxf, pdf->mixture.pdf_count - 1);

		out = pdf_generate(pdf->mixture.pdfs[idx]);
		break;
	}
	default:
		assert(0);
	}

	return out;
}

bool_t lambertian_scatter(const material_t *mtl, const ray_t *in_ray, const ray_hit_t *hit, vec3_t *attenuation,
                          pdf_t *pdf) {
	(void)in_ray;
	const vec3_t target = hit->normal + random_in_unit_sphere();
	*attenuation = mtl->albedo_texture.callback(&mtl->albedo_texture, hit->uv, hit->point);
	*pdf = pdf_init_cosine(hit->normal);
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
