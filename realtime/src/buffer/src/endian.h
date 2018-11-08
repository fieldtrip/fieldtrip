/*
 * Copyright (C) 2017, Robert Oostenveld
 * Donders Institute for Brain, Cognition and Behaviour; Radboud University; NL
 */

#ifndef ENDIAN_H_
#define ENDIAN_H_

#include <stdint.h>

#include "platform_includes.h"
#include "buffer.h"
#include "message.h"

#ifdef __cplusplus
extern "C" {
#endif

/* definition of functions to assist with big/little endian conversion, see endian.c */
void endian_swap16(unsigned int numel, void *data);
void endian_swap32(unsigned int numel, void *data);
void endian_swap64(unsigned int numel, void *data);
int endian_swap_buf_to_native(UINT16_T command, UINT32_T bufsize, void *buf);
int endian_swap_from_native(UINT16_T orgCommand, message_t *msg);

#ifdef __cplusplus
}
#endif

#endif
