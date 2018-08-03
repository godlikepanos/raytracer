#include <engine/misc/public/tga.hpp>

int main(int, char**)
{
	using namespace RT_NAMESPACE;

	const u32 width = 256;
	const u32 height = 128;
	u8 data[height][width][3];
	for(u32 w = 0; w < width; ++w)
		for(u32 h = 0; h < height; ++h)
		{
			u8 r = ((f32)w / width) * 255;
			u8 g = ((f32)h / height) * 255;
			u8 b = 0.2 * 255;

			data[h][w][0] = r;
			data[h][w][1] = g;
			data[h][w][2] = b;
		}

	save_tga(L"./image.tga", &data[0][0][0], sizeof(data), width, height);

	return 0;
}