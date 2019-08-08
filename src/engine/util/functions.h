#pragma once

#include <engine/util/std.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include <x86intrin.h>

// Mersenne Twister
void seed_mt(u32_t seed);

// Mersenne Twister
u32_t rand_mt(void);

static inline f32_t rand_0f_to_1f() {
#if 1
	const f32_t out = (f32_t)rand_mt() / 0xFFFFFFFF;
#elif 0
	const f32_t out = (f32_t)rand() / RAND_MAX;
#else
	unsigned rand_val;
	_rdrand32_step(&rand_val);
	const f32_t out = (f32_t)rand_val / 0xFFFFFFFF;
#endif
	return out;
}

#define MIN(x_, y_) ((x_) < (y_) ? (x_) : (y_))
#define MAX(x_, y_) ((x_) > (y_) ? (x_) : (y_))
