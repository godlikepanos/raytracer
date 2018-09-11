#include <engine/misc/public/tga.hpp>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>

RT_BEGIN_NAMESPACE

error_t save_tga(const char_t* fname, const u8_t* rgb, size_t rgb_size, u32_t width, u32_t height)
{
	assert(fname && rgb && rgb_size > 0 && (rgb_size % 3) == 0 && width && height);

	// Convert to utf8
	std::vector<char> ansi_string;
	ansi_string.resize(std::wcslen(fname) + 1);

	mbstate_t mbs;
	memset(&mbs, 0, sizeof(mbs));

	int ret = std::wcsrtombs(&ansi_string[0], &fname, ansi_string.size(), &mbs);
	assert(ret < (int)ansi_string.size());
	ansi_string[(u32_t)ret] = '\0';

	// Open file
	std::ofstream file;
	file.open(&ansi_string[0], std::ios::binary);

	if(!file.good())
		return ERROR_FILE_ACCESS;

	// Write header
	u8_t header[18] = {0};
	header[2] = 2;
	header[12] = width & 0xFF;
	header[13] = (width >> 8) & 0xFF;
	header[14] = height & 0xFF;
	header[15] = (height >> 8) & 0xFF;
	header[16] = 24;

	file.write((const char*)header, sizeof(header));

	// Write data
	for(u32_t y = 0; y < height; y++)
		for(u32_t x = 0; x < width; x++)
		{
			const u32_t idx = (u32_t)y * width * 3 + x * 3;
			assert(idx + 2 < rgb_size);

			file.put((char)rgb[idx + 2]);
			file.put((char)rgb[idx + 1]);
			file.put((char)rgb[idx + 0]);
		}

	// Write footer
	static const char footer[26] = "\0\0\0\0"
								   "\0\0\0\0"
								   "TRUEVISION-XFILE"
								   ".";
	file.write(footer, sizeof(footer));

	return ERROR_NONE;
}

RT_END_NAMESPACE
