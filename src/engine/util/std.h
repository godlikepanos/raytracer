#pragma once

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <wchar.h>

typedef uint8_t u8_t;
typedef uint16_t u16_t;
typedef uint32_t u32_t;
typedef uint64_t u64_t;
typedef int8_t i8_t;
typedef int16_t i16_t;
typedef int32_t i32_t;
typedef int64_t i64_t;
typedef float f32_t;
typedef double f64_t;
typedef char char_t;
typedef _Bool bool_t;
typedef int error_t;

static const error_t ERROR_NONE = 0;
static const error_t ERROR_FILE_ACCESS = 1;

static const bool_t TRUE = 1;
static const bool_t FALSE = 0;

#define NELEMS(arr) (sizeof(arr) / sizeof(arr[0]))
