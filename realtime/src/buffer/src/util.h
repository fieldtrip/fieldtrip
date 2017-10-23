/*
 * Copyright (C) 2017, Robert Oostenveld
 * Donders Institute for Brain, Cognition and Behaviour; Radboud University; NL
 */

#ifdef __cplusplus
extern "C" {
#endif

/* definition of various utility functions, see util.c */
unsigned int bufread(int s, void *buf, unsigned int numel);
unsigned int bufwrite(int s, const void *buf, unsigned int numel);
unsigned int append(void **buf1, unsigned int bufsize1, void *buf2, unsigned int bufsize2);
void check_datatypes();
unsigned int wordsize_from_type(UINT32_T data_type);
const ft_chunk_t *find_chunk(const void *buf, unsigned int offset0, unsigned int size, UINT32_T chunk_type);
int check_event_array(unsigned int size, const void *buf);

#ifdef __cplusplus
}
#endif

