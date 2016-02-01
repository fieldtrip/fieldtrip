// CalcMD5.c
// 128 bit MD5 checksum: file, string, byte stream
// This function calculates a 128 bit checksum for arrays or files.
// Digest = CalcMD5(Data, [InClass], [OutClass])
// INPUT:
//   Data:   Data array or file name. Either numerical or CHAR array.
//           Currently only files and arrays with up to 2^32 bytes (2.1GB) are
//           accepted.
//   InClass: String to declare the type of the 1st input.
//           Optional. Default: 'Char'.
//           'File': [Data] is a file name as string. The digest is calculated
//                   for this file.
//           'Char': [Data] is a char array to calculate the digest for. Only the
//                   ASCII part of the Matlab CHARs is used, such that the digest
//                   is the same as if the Matlab string is written to a file as
//                   UCHAR, e.g. with FWRITE.
//           'Unicode': All bytes of the input [Data] are used to calculate the
//                   digest. If [Data] has a numerical type, this method is
//                   applied ever.
//   OutClass: String, format of the output. Just the first character matters.
//           Optional, default: 'hex'.
//           'hex': [1 x 32] string as lowercase hexadecimal number.
//           'HEX': [1 x 32] string as lowercase hexadecimal number.
//           'Dec': [1 x 16] double vector with UINT8 values.
//           'Base64': [1 x 22] string, encoded to base 64 (A:Z,a:z,0:9,+,/).
//
// OUTPUT:
//   Digest: A 128 bit number is replied in a format depending on [OutClass].
//
// EXAMPLES:
//   Three methods to get the MD5 of a file:
//   1. Direct file access (recommended):
//     MD5 = CalcMD5(which('CalcMD5.m'), 'File')
//   2. Import the file to a CHAR array (binary mode for exact line breaks!):
//     FID = fopen(which('CalcMD5.m'), 'rb');
//     S   = fread(FID, inf, 'uchar=>char');
//     fclose(FID);
//     MD5 = CalcMD5(S, 'char')
//   3. Import file as a byte stream:
//     FID = fopen(which('CalcMD5.m'), 'rb');
//     S   = fread(FID, inf, 'uint8=>uint8');
//     fclose(FID);
//     MD5 = CalcMD5(S, 'unicode');  // 'unicode' can be omitted here
//
//   Test data:
//     CalcMD5(char(0:511), 'char', 'HEX')
//       => F5C8E3C31C044BAE0E65569560B54332
//     CalcMD5(char(0:511), 'unicode')
//       => 3484769D4F7EBB88BBE942BB924834CD
//
// Compile with:
//   mex -O CalcMD5.c
// On Linux the C99 comments must be considered (thanks Sebastiaan Breedveld):
//   mex -O CFLAGS="\$CFLAGS -std=C99" CalcMD5.c
//
// Tested: Matlab 6.5, 7.7, 7.8, WinXP, [UnitTest]
//         Compiler: BCC5.5, LCC2.4/3.8, OpenWatcom 1.8
// Author: Jan Simon, Heidelberg, (C) 2006-2010 J@n-Simon.De
// License: BSD. This program is based on:
//          RFC 1321, MD5 Message-Digest Algorithm, April 1992
//          RSA Data Security, Inc. MD5 Message Digest Algorithm
//          Modifications:
//          - Acceleration:  Unrolled loops. Compacted macros FF, GG, HH, II.
//          - Mex-interface: Input and output from and to Matlab.
//
// See also: CalcCRC32.
//
// Michael Kleder has published a Java call to compute the MD5 and SHA sums:
//   http://www.mathworks.com/matlabcentral/fileexchange/8944

/**********************************************************************
 ** Copyright (C) 1990, RSA Data Security, Inc. All rights reserved. **
 **                                                                  **
 ** License to copy and use this software is granted provided that   **
 ** it is identified as the "RSA Data Security, Inc. MD5 Message     **
 ** Digest Algorithm" in all material mentioning or referencing this **
 ** software or this function.                                       **
 **                                                                  **
 ** License is also granted to make and use derivative works         **
 ** provided that such works are identified as "derived from the RSA **
 ** Data Security, Inc. MD5 Message Digest Algorithm" in all         **
 ** material mentioning or referencing the derived work.             **
 **                                                                  **
 ** RSA Data Security, Inc. makes no representations concerning      **
 ** either the merchantability of this software or the suitability   **
 ** of this software for any particular purpose.  It is provided "as **
 ** is" without express or implied warranty of any kind.             **
 **                                                                  **
 ** These notices must be retained in any copies of any part of this **
 ** documentation and/or software.                                   **
 **********************************************************************
 */
 
/*
% $JRev: R5.00z V:025 Sum:/kHGslMmCpAS Date:17-Dec-2009 12:46:26 $
% $File: CalcMD5\CalcMD5.c $
% History:
% 011: 20-Oct-2006 20:50, [16 x 1] -> [1 x 16] replied as double.
% 012: 01-Nov-2006 23:10, BUGFIX: hex output for 'Hex' input now.
% 015: 02-Oct-2008 14:47, Base64 output.
% 017: 19-Oct-2008 22:33, Accept numerical arrays as byte stream.
% 023: 15-Dec-2009 16:53, BUGFIX: UINT32 has 32 bits on 64 bit systems now.
%      Thanks to Sebastiaan Breedveld!
*/

// Headers:
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "mex.h"

// Assume 32 bit array dimensions for Matlab 6.5:
// See option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef mwSize
#define mwSize  int
#define mwIndex int
#endif

// Types:
typedef unsigned char UCHAR;
typedef unsigned int  UINT;
typedef unsigned char * POINTER;   // generic pointer
typedef UINT32_T UINT32;           // four byte word (defined in tmwtypes.h)

typedef struct {
  UINT32 state[4];   // state (ABCD)
  UINT32 count[2];   // number of bits, modulo 64 (lsb first)
  UCHAR buffer[64];  // input buffer
} MD5_CTX;

// Prototypes:
void MD5Init     (MD5_CTX *);
void MD5Update   (MD5_CTX *, UCHAR *, UINT);
void MD5Final    (UCHAR[16], MD5_CTX *);
void MD5Transform(UINT32[4], UCHAR[64]);
void MD5Encode   (UCHAR *, UINT32 *, UINT);
void MD5Array    (UCHAR *data, mwSize N, UCHAR digest[16]);
void MD5File     (char *FileName, UCHAR digest[16]);
void MD5Char     (mxChar *data, mwSize N, UCHAR digest[16]);
void ToHex       (const UCHAR In[16], char *Out, int LowerCase);
void ToBase64    (const UCHAR In[16], char *Out);

// Constants for MD5Transform routine:
#define S11 7
#define S12 12
#define S13 17
#define S14 22
#define S21 5
#define S22 9
#define S23 14
#define S24 20
#define S31 4
#define S32 11
#define S33 16
#define S34 23
#define S41 6
#define S42 10
#define S43 15
#define S44 21

static UCHAR PADDING[64] = {
  0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

// F, G, H and I are basic MD5 functions:
#define F(x, y, z) (((x) & (y)) | ((~x) & (z)))
#define G(x, y, z) (((x) & (z)) | ((y) & (~z)))
#define H(x, y, z) ((x) ^ (y) ^ (z))
#define I(x, y, z) ((y) ^ ((x) | (~z)))

// ROTATE_LEFT rotates x left n bits:
// Rotation is separate from addition to prevent recomputation.
#define ROTATE_LEFT(x, n) (((x) << (n)) | ((x) >> (32 - (n))))

// FF, GG, HH, and II transformations for rounds 1, 2, 3, and 4:
#define FF(a, b, c, d, x, s, ac) { \
 (a) = ROTATE_LEFT((a) + F((b), (c), (d)) + (x) + (UINT32)(ac), (s)) + (b); }
#define GG(a, b, c, d, x, s, ac) { \
 (a) = ROTATE_LEFT((a) + G((b), (c), (d)) + (x) + (UINT32)(ac), (s)) + (b); }
#define HH(a, b, c, d, x, s, ac) { \
 (a) = ROTATE_LEFT((a) + H((b), (c), (d)) + (x) + (UINT32)(ac), (s)) + (b); }
#define II(a, b, c, d, x, s, ac) { \
 (a) = ROTATE_LEFT((a) + I((b), (c), (d)) + (x) + (UINT32)(ac), (s)) + (b); }

// Length of the file buffer (must be < 2^31 for INT conversion):
#define BUFFER_LEN 1024
static UCHAR buffer[BUFFER_LEN];

// MD5 initialization. Begins an MD5 operation, writing a new context. =========
void MD5Init(MD5_CTX *context)
{
  // Load magic initialization constants:
  context->count[0] = 0;
  context->count[1] = 0;
  context->state[0] = 0x67452301;
  context->state[1] = 0xefcdab89;
  context->state[2] = 0x98badcfe;
  context->state[3] = 0x10325476;
}

// MD5 block update operation. Continues an MD5 message-digest operation,
// processing another message block, and updating the context.
void MD5Update(MD5_CTX *context, UCHAR *input, UINT inputLen)
{
  UINT index, partLen;
  int i, inputLenM63;

  // Compute number of bytes mod 64:
  index = (UINT)((context->count[0] >> 3) & 0x3F);
  
  // Update number of bits:
  if ((context->count[0] += ((UINT32)inputLen << 3)) < ((UINT32)inputLen << 3)) {
    context->count[1]++;
  }
  context->count[1] += ((UINT32)inputLen >> 29);
  
  partLen = 64 - index;
  
  // Transform as many times as possible:
  if (inputLen >= partLen) {
    memcpy((POINTER)&context->buffer[index], (POINTER)input, partLen);
    MD5Transform(context->state, context->buffer);
    
    inputLenM63 = inputLen - 63;
    for (i = partLen; i < inputLenM63; i += 64) {
      MD5Transform(context->state, &input[i]);
    }
    
    // Buffer remaining input: index = 0
    memcpy((POINTER)&context->buffer[0], (POINTER)&input[i], inputLen - i);
  } else {
    // Buffer remaining input: i = 0
    memcpy((POINTER)&context->buffer[index], (POINTER)input, inputLen);
  }
  
  return;
}

// Finalize MD5: ===============================================================
// Ends an MD5 message-digest operation, writing the message digest and zeroing
// the context.
void MD5Final(UCHAR digest[16], MD5_CTX *context)
{
  UCHAR bits[8];
  UINT index, padLen;

  // Save number of bits:
  MD5Encode(bits, context->count, 2);

  // Pad out to 56 mod 64:
  index  = (UINT)((context->count[0] >> 3) & 0x3f);
  padLen = (index < 56) ? (56 - index) : (120 - index);
  MD5Update(context, PADDING, padLen);
  
  // Append length before padding:
  MD5Update(context, bits, 8);
  
  // Store state in digest:
  MD5Encode(digest, context->state, 4);
  
  // Zero sensitive information:
  memset((POINTER)context, 0, sizeof(MD5_CTX));
}

// MD5 basic transformation. Transforms state based on block: ==================
void MD5Transform(UINT32 state[4], UCHAR block[64])
{
  UINT32 a = state[0],
         b = state[1],
         c = state[2],
         d = state[3],
         x[16];

  // Unroll the loop for speed:
  // UINT i, j;
  // for (i = 0, j = 0; j < 64; i++, j += 4) {
  //   x[i] = ((UINT32)block[j]) | (((UINT32)block[j + 1]) << 8) |
  //          (((UINT32)block[j + 2]) << 16) | (((UINT32)block[j + 3]) << 24);
  // }
  x[0]  = ( (UINT32)block[0])         | (((UINT32)block[1])  << 8) |
          (((UINT32)block[2]) << 16)  | (((UINT32)block[3])  << 24);
  x[1]  = ( (UINT32)block[4])         | (((UINT32)block[5])  << 8) |
          (((UINT32)block[6]) << 16)  | (((UINT32)block[7])  << 24);
  x[2]  = ( (UINT32)block[8])         | (((UINT32)block[9])  << 8) |
          (((UINT32)block[10]) << 16) | (((UINT32)block[11]) << 24);
  x[3]  = ( (UINT32)block[12])        | (((UINT32)block[13]) << 8) |
          (((UINT32)block[14]) << 16) | (((UINT32)block[15]) << 24);
  x[4]  = ( (UINT32)block[16])        | (((UINT32)block[17]) << 8) |
          (((UINT32)block[18]) << 16) | (((UINT32)block[19]) << 24);
  x[5]  = ( (UINT32)block[20])        | (((UINT32)block[21]) << 8) |
          (((UINT32)block[22]) << 16) | (((UINT32)block[23]) << 24);
  x[6]  = ( (UINT32)block[24])        | (((UINT32)block[25]) << 8) |
          (((UINT32)block[26]) << 16) | (((UINT32)block[27]) << 24);
  x[7]  = ( (UINT32)block[28])        | (((UINT32)block[29]) << 8) |
          (((UINT32)block[30]) << 16) | (((UINT32)block[31]) << 24);
  x[8]  = ( (UINT32)block[32])        | (((UINT32)block[33]) << 8) |
          (((UINT32)block[34]) << 16) | (((UINT32)block[35]) << 24);
  x[9]  = ( (UINT32)block[36])        | (((UINT32)block[37]) << 8) |
          (((UINT32)block[38]) << 16) | (((UINT32)block[39]) << 24);
  x[10] = ( (UINT32)block[40])        | (((UINT32)block[41]) << 8) |
          (((UINT32)block[42]) << 16) | (((UINT32)block[43]) << 24);
  x[11] = ( (UINT32)block[44])        | (((UINT32)block[45]) << 8) |
          (((UINT32)block[46]) << 16) | (((UINT32)block[47]) << 24);
  x[12] = ( (UINT32)block[48])        | (((UINT32)block[49]) << 8) |
          (((UINT32)block[50]) << 16) | (((UINT32)block[51]) << 24);
  x[13] = ( (UINT32)block[52])        | (((UINT32)block[53]) << 8) |
          (((UINT32)block[54]) << 16) | (((UINT32)block[55]) << 24);
  x[14] = ( (UINT32)block[56])        | (((UINT32)block[57]) << 8) |
          (((UINT32)block[58]) << 16) | (((UINT32)block[59]) << 24);
  x[15] = ( (UINT32)block[60])        | (((UINT32)block[61]) << 8) |
          (((UINT32)block[62]) << 16) | (((UINT32)block[63]) << 24);
  
  // Round 1
  FF(a, b, c, d, x[ 0], S11, 0xd76aa478);  // 1
  FF(d, a, b, c, x[ 1], S12, 0xe8c7b756);  // 2
  FF(c, d, a, b, x[ 2], S13, 0x242070db);  // 3
  FF(b, c, d, a, x[ 3], S14, 0xc1bdceee);  // 4
  FF(a, b, c, d, x[ 4], S11, 0xf57c0faf);  // 5
  FF(d, a, b, c, x[ 5], S12, 0x4787c62a);  // 6
  FF(c, d, a, b, x[ 6], S13, 0xa8304613);  // 7
  FF(b, c, d, a, x[ 7], S14, 0xfd469501);  // 8
  FF(a, b, c, d, x[ 8], S11, 0x698098d8);  // 9
  FF(d, a, b, c, x[ 9], S12, 0x8b44f7af);  // 10
  FF(c, d, a, b, x[10], S13, 0xffff5bb1);  // 11
  FF(b, c, d, a, x[11], S14, 0x895cd7be);  // 12
  FF(a, b, c, d, x[12], S11, 0x6b901122);  // 13
  FF(d, a, b, c, x[13], S12, 0xfd987193);  // 14
  FF(c, d, a, b, x[14], S13, 0xa679438e);  // 15
  FF(b, c, d, a, x[15], S14, 0x49b40821);  // 16

  // Round 2
  GG(a, b, c, d, x[ 1], S21, 0xf61e2562);  // 17
  GG(d, a, b, c, x[ 6], S22, 0xc040b340);  // 18
  GG(c, d, a, b, x[11], S23, 0x265e5a51);  // 19
  GG(b, c, d, a, x[ 0], S24, 0xe9b6c7aa);  // 20
  GG(a, b, c, d, x[ 5], S21, 0xd62f105d);  // 21
  GG(d, a, b, c, x[10], S22,  0x2441453);  // 22
  GG(c, d, a, b, x[15], S23, 0xd8a1e681);  // 23
  GG(b, c, d, a, x[ 4], S24, 0xe7d3fbc8);  // 24
  GG(a, b, c, d, x[ 9], S21, 0x21e1cde6);  // 25
  GG(d, a, b, c, x[14], S22, 0xc33707d6);  // 26
  GG(c, d, a, b, x[ 3], S23, 0xf4d50d87);  // 27

  GG(b, c, d, a, x[ 8], S24, 0x455a14ed);  // 28
  GG(a, b, c, d, x[13], S21, 0xa9e3e905);  // 29
  GG(d, a, b, c, x[ 2], S22, 0xfcefa3f8);  // 30
  GG(c, d, a, b, x[ 7], S23, 0x676f02d9);  // 31
  GG(b, c, d, a, x[12], S24, 0x8d2a4c8a);  // 32

  // Round 3
  HH(a, b, c, d, x[ 5], S31, 0xfffa3942);  // 33
  HH(d, a, b, c, x[ 8], S32, 0x8771f681);  // 34
  HH(c, d, a, b, x[11], S33, 0x6d9d6122);  // 35
  HH(b, c, d, a, x[14], S34, 0xfde5380c);  // 36
  HH(a, b, c, d, x[ 1], S31, 0xa4beea44);  // 37
  HH(d, a, b, c, x[ 4], S32, 0x4bdecfa9);  // 38
  HH(c, d, a, b, x[ 7], S33, 0xf6bb4b60);  // 39
  HH(b, c, d, a, x[10], S34, 0xbebfbc70);  // 40
  HH(a, b, c, d, x[13], S31, 0x289b7ec6);  // 41
  HH(d, a, b, c, x[ 0], S32, 0xeaa127fa);  // 42
  HH(c, d, a, b, x[ 3], S33, 0xd4ef3085);  // 43
  HH(b, c, d, a, x[ 6], S34,  0x4881d05);  // 44
  HH(a, b, c, d, x[ 9], S31, 0xd9d4d039);  // 45
  HH(d, a, b, c, x[12], S32, 0xe6db99e5);  // 46
  HH(c, d, a, b, x[15], S33, 0x1fa27cf8);  // 47
  HH(b, c, d, a, x[ 2], S34, 0xc4ac5665);  // 48

  // Round 4
  II(a, b, c, d, x[ 0], S41, 0xf4292244);  // 49
  II(d, a, b, c, x[ 7], S42, 0x432aff97);  // 50
  II(c, d, a, b, x[14], S43, 0xab9423a7);  // 51
  II(b, c, d, a, x[ 5], S44, 0xfc93a039);  // 52
  II(a, b, c, d, x[12], S41, 0x655b59c3);  // 53
  II(d, a, b, c, x[ 3], S42, 0x8f0ccc92);  // 54
  II(c, d, a, b, x[10], S43, 0xffeff47d);  // 55
  II(b, c, d, a, x[ 1], S44, 0x85845dd1);  // 56
  II(a, b, c, d, x[ 8], S41, 0x6fa87e4f);  // 57
  II(d, a, b, c, x[15], S42, 0xfe2ce6e0);  // 58
  II(c, d, a, b, x[ 6], S43, 0xa3014314);  // 59
  II(b, c, d, a, x[13], S44, 0x4e0811a1);  // 60
  II(a, b, c, d, x[ 4], S41, 0xf7537e82);  // 61
  II(d, a, b, c, x[11], S42, 0xbd3af235);  // 62
  II(c, d, a, b, x[ 2], S43, 0x2ad7d2bb);  // 63
  II(b, c, d, a, x[ 9], S44, 0xeb86d391);  // 64

  state[0] += a;
  state[1] += b;
  state[2] += c;
  state[3] += d;

  memset((POINTER)x, 0, sizeof(x));
}

// Encodes input (UINT32) into output (UCHAR) (length is divided by 4) =========
void MD5Encode(UCHAR *output, UINT32 *input, UINT len)
{
  UINT j;
  
  for (j = 0; j < len; j++) {
    *output++ = (UCHAR)( *input          & 0xff);
    *output++ = (UCHAR)((*input   >>  8) & 0xff);
    *output++ = (UCHAR)((*input   >> 16) & 0xff);
    *output++ = (UCHAR)((*input++ >> 24) & 0xff);
  }
}

// Calcualte digest: ===========================================================
void MD5Char(mxChar *array, mwSize inputLen, UCHAR digest[16])
{
  // Process string: Matlab stores strings as mxChar, which are 2 bytes per
  // character. This function considers the first byte of each CHAR only, which
  // is equivalent to calculate the sum after a conversion to a ASCII UCHAR
  // string.
  MD5_CTX context;
  UINT Chunk;
  UCHAR *bufferP, *bufferEnd = buffer + BUFFER_LEN, *arrayP;
  
  // Limit length to 32 bit address, because I cannot test this function
  // with 64 bit arrays currently (under construction):
  if (inputLen >> 31 != 0) {  // Detect sign-bit if mwSize is int
     mexErrMsgTxt("*** CalcMD5[mex]: Input > 2^31 byte not handled yet.");
  }
  
  arrayP = (UCHAR *) array;  // UCHAR *, not mxChar *!
  
  MD5Init(&context);
  
  // Copy chunks of input data - only the first byte of each mxChar:
  Chunk = inputLen / BUFFER_LEN;
  while (Chunk--) {
     bufferP = buffer;
     while (bufferP < bufferEnd) {
        *bufferP++ = *arrayP;
        arrayP    += 2;
     }
     
     MD5Update(&context, buffer, BUFFER_LEN);
  }
  
  // Last chunk:
  Chunk = inputLen % BUFFER_LEN;
  if (Chunk != 0) {
     bufferEnd = buffer + Chunk;
     bufferP   = buffer;
     while (bufferP < bufferEnd) {
        *bufferP++ = *arrayP;
        arrayP    += 2;
     }
     
     MD5Update(&context, buffer, Chunk);
  }
  
  MD5Final(digest, &context);
  
  return;
}

// Array of any type as byte stream: ===========================================
void MD5Array(UCHAR *array, mwSize inputLen, UCHAR digest[16])
{
  MD5_CTX context;
  
  // Limit length to 32 bit address, because I cannot test this function
  // with 64 bit arrays currently (under construction):
  if (inputLen >> 31 != 0) {  // Detect sign-bit if mwSize is signed int
     mexErrMsgTxt("*** CalcMD5[mex]: Input > 2^31 byte not handled yet.");
  }
  
  MD5Init(&context);
  MD5Update(&context, array, (UINT) inputLen);
  MD5Final(digest, &context);
}

// File as byte stream: ========================================================
void MD5File(char *filename, UCHAR digest[16])
{
  FILE *FID;
  MD5_CTX context;
  int len;
  UINT32 allLen = 0;
  
  // Open the file in binary mode:
  if ((FID = fopen(filename, "rb")) == NULL) {
     mexPrintf("*** Error for file: [%s]\n", filename);
     mexErrMsgTxt("*** CalcMD5[mex]: Cannot open file.");
  }
  
  MD5Init(&context);
  while ((len = fread(buffer, 1, BUFFER_LEN, FID)) != 0) {
     // Limit length to 32 bit address, because I cannot test this function
     // with 64 bit arrays currently (under construction):
     allLen += len;
     if (allLen > 2147483647) {  // 2^31
        fclose(FID);
        mexErrMsgTxt("*** CalcMD5[mex]: Cannot handle files > 2.1GB yet.");
     }
     
     MD5Update(&context, buffer, (UINT) len);
  }
  MD5Final(digest, &context);

  fclose(FID);
}

// Output of 16 UCHARs as 32 character hexadecimals: ===========================
void ToHex(const UCHAR digest[16], char *output, int LowerCase)
{
  char *outputEnd;
  
  if (LowerCase) {
    for (outputEnd = output + 32; output < outputEnd; output += 2) {
      sprintf(output, "%02x", *(digest++));
    }
  } else {  // Upper case:
    for (outputEnd = output + 32; output < outputEnd; output += 2) {
      sprintf(output, "%02X", *(digest++));
    }
  }
  
  return;
}

// BASE64 encoded output: ======================================================
void ToBase64(const UCHAR In[16], char *Out)
{
   // The base64 encoded string is shorter than the hex string.
   // Needed length: ((len + 2) / 3 * 4) + 1, here fixed to 22+1 here (trailing
   // 0 included).
   static const UCHAR B64[] =
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

   int i;
   char *p;
   const UCHAR *s;
   
   p = Out;
   s = In;
   for (i = 0; i < 5; i++) {
      *p++ = B64[(*s >> 2) & 0x3F];
      *p++ = B64[((*s & 0x3) << 4)   | ((s[1] & 0xF0) >> 4)];
      *p++ = B64[((s[1] & 0xF) << 2) | ((s[2] & 0xC0) >> 6)];
      *p++ = B64[s[2] & 0x3F];
      s   += 3;
   }
   
   *p++ = B64[(*s >> 2) & 0x3F];
   *p++ = B64[((*s & 0x3) << 4)];
   *p   = '\0';
   
   return;
}

// Main function: ==============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Mex interface:
  // - Define default values of optional arguments.
  // - Forward input data to different calculators according to the input type.
  // - Convert digest to output format.
  
  char   *FileName, InType, hexOut[33], b64Out[23];
  UCHAR  digest[16], *digestP, OutType = 'h';
  int    isFile = false, isUnicode = false;
  double *outP, *outEnd;
  
  // Check number of inputs and outputs:
  if (nrhs == 0 || nrhs > 3) {
    mexErrMsgTxt("*** CalcMD5[mex]: 1 to 3 inputs required.");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("*** CalcMD5[mex]: Too many output arguments.");
  }
 
  // If 2nd input starts with 'f', treat string in 1st argument as file name:
  if (nrhs >= 2 && mxGetNumberOfElements(prhs[1]) > 0) {
    if (mxIsChar(prhs[1]) == 0) {
      mexErrMsgTxt("*** CalcMD5[mex]: 2nd input must be a string.");
    }
    
    InType    = (char) tolower(*(POINTER) mxGetData(prhs[1]));
    isFile    = (InType == 'f');
    isUnicode = (InType == 'u');
  }  // Default otherwise!
  
  // Output type - default: hex:
  if (nrhs == 3 && !mxIsEmpty(prhs[2])) {
    if (mxIsChar(prhs[2]) == 0) {
      mexErrMsgTxt("*** CalcMD5[mex]: 3rd input must be a string.");
    }
    
    OutType = *(POINTER) mxGetData(prhs[2]);  // Just 1st character
  }
     
  // Calculate check sum:
  if (isFile) {
     if ((FileName = mxArrayToString(prhs[0])) == NULL) {
        mexErrMsgTxt("*** CalcMD5[mex]: Cannot get file name.");
     }
     MD5File(FileName, digest);
     mxFree(FileName);
     
  } else if (mxIsNumeric(prhs[0]) || isUnicode) {
     MD5Array((POINTER) mxGetData(prhs[0]),
              mxGetNumberOfElements(prhs[0]) * mxGetElementSize(prhs[0]),
              digest);
              
  } else if (mxIsChar(prhs[0])) {
     MD5Char((mxChar *) mxGetData(prhs[0]),
             mxGetNumberOfElements(prhs[0]),
             digest);
     
  } else {
     mexErrMsgTxt("*** CalcMD5[mex]: Input type not accepted.");
  }
  
  // Create output:
  switch (OutType) {
    case 'H':
    case 'h':  // Hexadecimal upper/lower case:
      ToHex(digest, hexOut, OutType == 'h');
      plhs[0] = mxCreateString(hexOut);
      break;
      
    case 'D':
    case 'd':  // DOUBLE with integer values:
      plhs[0] = mxCreateDoubleMatrix(1, 16, mxREAL);
      outP    = mxGetPr(plhs[0]);
      digestP = digest;
      for (outEnd = outP + 16; outP < outEnd; outP++) {
        *outP = (double) *digestP++;
      }
      break;
    
    case 'B':
    case 'b':  // Base64:
      //strtobase64(b64Out, 26, digest, 16);  // included in LCC3.8
      //b64Out[24] = '\0';
      ToBase64(digest, b64Out);               // Locally implemented
      plhs[0] = mxCreateString(b64Out);
      break;
      
    default:
      mexErrMsgTxt("*** CalcMD5[mex]: Unknown output type.");
  }
  
  return;
}
