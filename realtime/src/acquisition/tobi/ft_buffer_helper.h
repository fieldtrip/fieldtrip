#ifndef FT_BUFFER_HELPER
#define FT_BUFFER_HELPER

#include "buffer.h"

// Forward declarations for high-level helpers functions
int ft_put_data(int ft_buffer_handle, int nchannels, int nsamples, const float
  *chan_samp);

int ft_put_hdr(int ft_buffer_handle, int nchann, int fsample, ft_chunk_t
  *chunk);

ft_chunk_t *ft_create_chanlab_chunk(int n, const char **labels);

void ft_buffer_serve(int port);

#endif
