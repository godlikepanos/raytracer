#pragma once

#include <engine/util/std.h>

error_t save_tga(const char_t* fname, const u8_t* rgb, size_t rgb_size, u32_t width, u32_t height);
