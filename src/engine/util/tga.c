#include <engine/util/tga.h>

error_t save_tga(const char_t *fname, const u8_t *rgb, size_t rgb_size, u32_t width, u32_t height) {
	assert(fname && rgb && rgb_size > 0 && (rgb_size % 3) == 0 && width && height);

	// Open file
	FILE *file = fopen(fname, "wb");

	if(!file)
		return ERROR_FILE_ACCESS;

	// Write header
	u8_t header[18] = {0};
	header[2]       = 2;
	header[12]      = width & 0xFF;
	header[13]      = (width >> 8) & 0xFF;
	header[14]      = height & 0xFF;
	header[15]      = (height >> 8) & 0xFF;
	header[16]      = 24;

	fwrite((const char *)header, 1, sizeof(header), file);

	// Write data
	for(u32_t y = 0; y < height; y++) {
		for(u32_t x = 0; x < width; x++) {
			const u32_t idx = (u32_t)y * width * 3 + x * 3;
			assert(idx + 2 < rgb_size);

			fputc((char)rgb[idx + 2], file);
			fputc((char)rgb[idx + 1], file);
			fputc((char)rgb[idx + 0], file);
		}
	}

	// Write footer
	static const char footer[26] = "\0\0\0\0"
	                               "\0\0\0\0"
	                               "TRUEVISION-XFILE"
	                               ".";
	fwrite(footer, 1, sizeof(footer), file);
	fclose(file);

	return ERROR_NONE;
}
