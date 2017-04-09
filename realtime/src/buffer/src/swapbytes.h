/*
 * Copyright (C) 2010, Stefan Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#ifndef SWAP_BYTES_H
#define SWAP_BYTES_H

#include <stdio.h>

#if defined(_MSC_VER)
#include <win32/stdint.h>
#elif defined(__BORLANDC__)
  /* without the following, compilation with the Borland command line tools fails -- SK */
  typedef          __int8     int8_t;
  typedef          __int16    int16_t;
  typedef          __int32    int32_t;
  typedef          __int64    int64_t;
  typedef unsigned __int8    uint8_t;
  typedef unsigned __int16   uint16_t;
  typedef unsigned __int32   uint32_t;
  typedef unsigned __int64   uint64_t;
#else
#include <stdint.h>
#endif

#define SwapTwoBytes(data)   ( (((data) >> 8) & 0x00FF) | (((data) << 8) & 0xFF00) )
#define SwapFourBytes(data)  ( (((data) >> 24) & 0x000000FF) | (((data) >> 8) & 0x0000FF00) | \ (((data) << 8) & 0x00FF0000) | (((data) << 24) & 0xFF000000) )
#define SwapEightBytes(data) ( (((data) >> 56) & 0x00000000000000FF) | (((data) >> 40) & 0x000000000000FF00) | \ (((data) >> 24) & 0x0000000000FF0000) | (((data) >> 8) & 0x00000000FF000000) | \ (((data) << 8) & 0x000000FF00000000) | (((data) << 24) & 0x0000FF0000000000) | \ (((data) << 40) & 0x00FF000000000000)

#if CPU_ARCHITECTURE == BIG_ENDIAN
#else
#endif

double   swap_double(double   x);
uint64_t swap_uint64(uint64_t x);
uint32_t swap_uint32(uint32_t x);
uint16_t swap_uint16(uint16_t x);

#endif
