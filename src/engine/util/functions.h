#pragma once

#include <engine/util/std.h>
#include <stdlib.h>
#include <time.h>

static inline f32_t rand_0f_to_1f() {
	return (f32_t)rand() / RAND_MAX;
}
