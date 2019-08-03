#include <engine/rt.h>

int main(int argc, char** argv) {
	const u32_t width  = 512;
	const u32_t height = 256;

	const vec3_t cam_pos	= {0.0f, 0.0f, 0.0f};
	const vec3_t eye		= {0.0f, 0.0f, -1.0f};
	const vec3_t up			= {0.0f, 1.0f, 0.0f};
	const mat4_t view_mat   = look_at(eye, cam_pos, up);
	const mat4_t cam_trf	= mat4_invert(&view_mat);
	const mat4_t proj_mat   = perspective(to_rad(60.0f), (f32_t)width / (f32_t)height, 0.1f, 100.0f);
	const mat4_t vp_mat		= mat4_mul_mat4(&proj_mat, &view_mat);
	const mat4_t inv_vp_mat = mat4_invert(&vp_mat);

	u8_t data[height][width][3];
	for(u32_t w = 0; w < width; ++w)
		for(u32_t h = 0; h < height; ++h) {
			vec2_t ndc = {(w + 0.5f) / width, (h + 0.5f) / height};
			ndc		   = ndc * vec2_init_f(2.0f) - vec2_init_f(1.0f);

			const vec3_t view_dir = mat4_mul_vec4(&inv_vp_mat, vec4_init_4f(ndc.x, ndc.y, 1.0f, 1.0f)).xyz;
			mat3_t		 cam_rot  = mat3_init_mat4(&cam_trf);
			ray_t		 ray;
			ray.direction = mat3_mul_vec3(&cam_rot, view_dir);
			ray.origin	= cam_pos;

			const vec3_t c_a   = {0.0f, 0.0f, 1.0f};
			const vec3_t c_b   = vec3_init_f(1.0f);
			vec3_t		 color = vec3_mix(c_a, c_b, (f32_t)h / height);

			ray_hit_t hit;
			sphere_t  s = sphere_init(vec3_init_3f(0.0f, 0.0f, -5.0f), 0.5f);
			if(ray_cast_sphere(&ray, &s, 0.01f, 100.f, &hit)) {
				color = hit.normal;
			}

			color		  = vec3_clamp(color, vec3_init_f(0.0f), vec3_init_f(1.0f));
			data[h][w][0] = color.x * 255;
			data[h][w][1] = color.y * 255;
			data[h][w][2] = color.z * 255;
		}

	save_tga("./image.tga", &data[0][0][0], sizeof(data), width, height);

	return 0;
}
