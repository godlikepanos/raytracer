#pragma once

#include <cstdint>
#include <cassert>

#define RT_NAMESPACE rt
#define RT_BEGIN_NAMESPACE \
	namespace RT_NAMESPACE \
	{
#define RT_END_NAMESPACE }

RT_BEGIN_NAMESPACE

using u8_t = uint8_t;
using u16_t = uint16_t;
using u32_t = uint32_t;
using u64_t = uint64_t;

using s8_t = int8_t;
using s16_t = int16_t;
using s32_t = int32_t;
using s64_t = int64_t;

using size_t = std::size_t;

using f32_t = float;
using f64_t = double;

using error_t = int;
const error_t ERROR_NONE = 0;
const error_t ERROR_FILE_ACCESS = 1;

using char_t = wchar_t;
using bool_t = bool;

RT_END_NAMESPACE
