#include "swapbytes.h"

/* Test platform Endianness */
int bigendian(void) {
	register int i;
	union {
		char Array[4];
		long Chars;
	} TestUnion;
	char c = 'a';
	for(i=0; i<4; i++)
		TestUnion.Array[i] = c++; /* FIXME: ++ on non-numeric type. */
	if (TestUnion.Chars == 0x61626364)
		return 1;
	else
		return 0;
}

double swap_double(double x) {
	double y;
    ((uint8_t*)(&y))[0] = ((uint8_t*)(&x))[7];	
    ((uint8_t*)(&y))[1] = ((uint8_t*)(&x))[6];	
    ((uint8_t*)(&y))[2] = ((uint8_t*)(&x))[5];	
    ((uint8_t*)(&y))[3] = ((uint8_t*)(&x))[4];	
    ((uint8_t*)(&y))[4] = ((uint8_t*)(&x))[3];	
    ((uint8_t*)(&y))[5] = ((uint8_t*)(&x))[2];	
    ((uint8_t*)(&y))[6] = ((uint8_t*)(&x))[1];	
    ((uint8_t*)(&y))[7] = ((uint8_t*)(&x))[0];	
	return y;
}

uint64_t swap_uint64(uint64_t x) {
	uint64_t y;
    ((uint8_t*)(&y))[0] = ((uint8_t*)(&x))[7];	
    ((uint8_t*)(&y))[1] = ((uint8_t*)(&x))[6];	
    ((uint8_t*)(&y))[2] = ((uint8_t*)(&x))[5];	
    ((uint8_t*)(&y))[3] = ((uint8_t*)(&x))[4];	
    ((uint8_t*)(&y))[4] = ((uint8_t*)(&x))[3];	
    ((uint8_t*)(&y))[5] = ((uint8_t*)(&x))[2];	
    ((uint8_t*)(&y))[6] = ((uint8_t*)(&x))[1];	
    ((uint8_t*)(&y))[7] = ((uint8_t*)(&x))[0];	
	return y;
}

uint32_t swap_uint32(uint32_t x) {
	uint32_t y;
    ((uint8_t*)(&y))[0] = ((uint8_t*)(&x))[3];	
    ((uint8_t*)(&y))[1] = ((uint8_t*)(&x))[2];	
    ((uint8_t*)(&y))[2] = ((uint8_t*)(&x))[1];	
    ((uint8_t*)(&y))[3] = ((uint8_t*)(&x))[0];	
	return y;
}

uint16_t swap_uint16(uint16_t x) {
	uint16_t y;
    ((uint8_t*)(&y))[0] = ((uint8_t*)(&x))[1];	
    ((uint8_t*)(&y))[1] = ((uint8_t*)(&x))[0];	
	return y;
}

