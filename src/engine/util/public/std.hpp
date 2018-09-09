#pragma once

#include <cstdint>
#include <cassert>

#define RT_NAMESPACE rt
#define RT_BEGIN_NAMESPACE \
	namespace RT_NAMESPACE \
	{
#define RT_END_NAMESPACE }

RT_BEGIN_NAMESPACE

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

using s8 = int8_t;
using s16 = int16_t;
using s32 = int32_t;
using s64 = int64_t;

using tsize = std::size_t;

using f32 = float;
using f64 = double;

using error = int;
const error ERROR_NONE = 0;
const error ERROR_FILE_ACCESS = 1;

using tchar = wchar_t;
using boolean = int;

RT_END_NAMESPACE
