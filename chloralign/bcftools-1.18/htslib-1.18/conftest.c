#ifdef __x86_64__
#include "x86intrin.h"
int main(int argc, char **argv) {
    __m128i a = _mm_set_epi32(1, 2, 3, 4), b = _mm_set_epi32(4, 3, 2, 1);
    __m128i c = _mm_shuffle_epi8(a, b);
    return *((char *) &c);
}
#else
int main(int argc, char **argv) { return 0; }
#endif
