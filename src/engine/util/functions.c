#include <engine/util/functions.h>

#define N (624)                               // length of state vector
#define M (397)                               // a period parameter
#define K (0x9908B0DF)                        // a magic constant
#define hibit(u) ((u)&0x80000000)             // mask all but highest   bit of u
#define lobit(u) ((u)&0x00000001)             // mask all but lowest    bit of u
#define lobits(u) ((u)&0x7FFFFFFF)            // mask     the highest   bit of u
#define mix_bits(u, v) (hibit(u) | lobits(v)) // move hi bit of u to hi bit of v

static __thread u32_t state[N + 1]; // state vector + 1 extra to not violate ANSI C
static __thread u32_t *next;        // next random value is computed from here
static __thread i32_t left = -1;    // can *next++ this many times before reloading

void seed_mt(u32_t seed) {

	u32_t x = (seed | 1U) & 0xFFFFFFFF;
	u32_t *s = state;
	i32_t j;

	for(left = 0, *s++ = x, j = N; --j; *s++ = (x *= 69069) & 0xFFFFFFFF)
		;
}

static u32_t reload_mt(void) {
	u32_t *p0 = state;
	u32_t *p2 = state + 2;
	u32_t *pM = state + M;
	u32_t s0;
	u32_t s1;
	i32_t j;

	if(left < -1) {
		seed_mt(4357);
	}

	left = N - 1, next = state + 1;

	for(s0 = state[0], s1 = state[1], j = N - M + 1; --j; s0 = s1, s1 = *p2++) {
		*p0++ = *pM++ ^ (mix_bits(s0, s1) >> 1) ^ (lobit(s1) ? K : 0);
	}

	for(pM = state, j = M; --j; s0 = s1, s1 = *p2++) {
		*p0++ = *pM++ ^ (mix_bits(s0, s1) >> 1) ^ (lobit(s1) ? K : 0);
	}

	s1 = state[0], *p0 = *pM ^ (mix_bits(s0, s1) >> 1) ^ (lobit(s1) ? K : 0);
	s1 ^= (s1 >> 11);
	s1 ^= (s1 << 7) & 0x9D2C5680;
	s1 ^= (s1 << 15) & 0xEFC60000;
	return (s1 ^ (s1 >> 18));
}

u32_t rand_mt(void) {
	u32_t y;

	if(--left < 0) {
		return (reload_mt());
	}

	y = *next++;
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9D2C5680U;
	y ^= (y << 15) & 0xEFC60000U;
	return y ^ (y >> 18);
}
