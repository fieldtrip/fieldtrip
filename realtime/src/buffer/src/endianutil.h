/*
 * Copyright (C) 2017, Robert Oostenveld
 * Donders Institute for Brain, Cognition and Behaviour; Radboud University; NL
 */

#ifdef __cplusplus
extern "C" {
#endif

/* definition of functions to assist with big/little endian conversion, see endianutil.c */
void ft_swap16(unsigned int numel, void *data);
void ft_swap32(unsigned int numel, void *data);
void ft_swap64(unsigned int numel, void *data);
int ft_swap_buf_to_native(UINT16_T command, UINT32_T bufsize, void *buf);
int ft_convert_chunks_from_native(UINT32_T size, UINT32_T nchans, void *buf);
int ft_swap_from_native(UINT16_T orgCommand, message_t *msg);

#ifdef __cplusplus
}
#endif

