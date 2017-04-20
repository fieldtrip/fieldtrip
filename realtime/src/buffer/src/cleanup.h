/*
 * Copyright (C) 2017, Robert Oostenveld
 * Donders Institute for Brain, Cognition and Behaviour; Radboud University; NL
 */

#ifdef __cplusplus
extern "C" {
#endif

/* definition of the functions used for thread cancelation, see cleanup.c */
void cleanup_message(void **arg);
void cleanup_header(void **arg);
void cleanup_data(void **arg);
void cleanup_event(void **arg);
void cleanup_buf(void **arg);
void cleanup_socket(int *);

#ifdef __cplusplus
}
#endif

