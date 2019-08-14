#pragma once

#include <unistd.h>

static inline u32_t get_gpu_core_count() {
	return sysconf(_SC_NPROCESSORS_ONLN);
}
