#pragma once

#include <engine/util/public/std.hpp>

RT_BEGIN_NAMESPACE

error save_tga(const tchar* fname, const u8* rgb, tsize rgb_size, u32 width, u32 height);

RT_END_NAMESPACE
