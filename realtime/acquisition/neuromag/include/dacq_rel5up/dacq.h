/* Automatically generated file, chages will not survive!*/
/* Copyright (c) 2009 Elekta Neuromag Oy, all rights reserved.*/
/* Generated  Tue Dec 29 12:48:12 EET 2009  by  skesti @ tammi.neuromag.fi  */

#ifndef dacq_h_included
#define dacq_h_included
/*
 *
 * Copyright 1990
 *
 * Matti Hamalainen
 * Low Temperature Laboratory
 * Helsinki University of Technology
 * SF-02150 Espoo
 * FINLAND
 *
 * No part of this program may be photocopied, reproduced,
 * or translated to another program language without the
 * prior written consent of the author.
 *
 * $Header: /sw/tmp/build/libs-dacq-1.2.11/src/RCS/dacq.hh,v 1.4 2009/01/12 11:00:38 skesti Exp $
 * $Id$
 */

#include "dacq_clients.h"
#include "dacq_sem.h"
#include "dacq_shmem.h"
#include <fiff.h>
#include <stdarg.h>

#define DACQ_MAGN_CH(x) 10000 + ((x) % 10000)
#define DACQ_EL_CH(x)   20000 + ((x) % 10000)
#define DACQ_STIM_CH(x) 30000 + ((x) % 10000)

typedef unsigned short u_pack;

typedef void (* Dacq_async_finish)(int *sock, char *command, char *response);

#if defined(__cplusplus)
extern "C" {
#endif
/* client_socket.c */
extern int dacq_connect_client(int id);
extern int dacq_disconnect_client(int *sock, int id);
extern int dacq_client_command(int *sock, int id, int cmd);
extern void dacq_set_data_filter(int *kinds, int nkind);
extern int dacq_client_receive_tag(int *sock, int id);
/* dummy_process_tag.c */
extern int (*dacq_client_process_tag)(fiffTag);
/* remote_clients.c */
extern void dacq_set_remote_timeouts(int new_read, int new_write);
extern int dacq_remote_async_find(int *sock);
extern void dacq_remote_async_add(int *sock, Dacq_async_finish finish);
extern void dacq_remote_async_remove(int *sock);
extern void dacq_remote_command_close(int *sock);
extern char *dacq_receive_reply_new(int *sock);
extern int dacq_connect_to_service(char *host, char *service);
extern int dacq_connect_collector(char *host);
extern int dacq_connect_isotrak(char *host);
extern int dacq_connect_xplotter(char *host);
extern char *dacq_remote_command_new(int *sock, char *command);
extern int dacq_remote_command_send_new(int *sock, char *command);
extern int dacq_collector_status_new(int *sock);
extern int dacq_remote_command_login(int *sock);
extern int dacq_remote_command_logof(int *sock);
/* xclient_socket.c */
extern int dacq_disconnect_xclient(int *sock, int id);
/* pack.c */
extern void dacq_pack_3852A(float *orig, u_pack *packed, int nrd);
extern void dacq_unpack_3852A(u_pack *packed, int nrd, float *orig);
extern void dacq_unpack_3852A_16(short *packed, float *ranges, int nchan, int nsamp, float *orig);
extern void dacq_unpack_3852A_16_t(short *packed, float *ranges, int nchan, int nsamp, float **orig);
extern void dacq_unpack_3852A_16_sel(short *packed, float *ranges, int nchan, int nsamp, int *selects, int nsel, float *orig);
extern void dacq_unpack_3852A_16_sel_t(short *packed, float *ranges, int nchan, int nsamp, int *selects, int nsel, float **orig);
extern void dacq_unpack_3852A_16_sel_ch_info_t(short *packed, fiffChInfo chs, int nchan, int nsamp, int *selects, int nsel, float **orig);
/* variables.c */
extern void dacq_set_constant_modify(int state);
extern void set_constant_modify(int state);
extern void dacq_set_var_prefix(char *new_prefix);
extern void dacq_set_ignore_prefix(char *new_prefix);
extern int dacq_set_int_var(char *name, int val);
extern int dacq_set_float_var(char *name, double val);
extern int dacq_set_text_var(char *name, char *val);
extern int dacq_append_text_var(char *name, char *val);
extern int dacq_get_int_var(char *name, int *val);
extern int dacq_get_float_var(char *name, float *val);
extern int dacq_get_text_var(char *name, char **val);
extern char *dacq_get_formatted_var(char *name);
extern void dacq_sort_variables(void);
extern int dacq_list_variables(FILE *out);
extern int dacq_define_variable(char *def);
extern int dacq_hide_variable(char *name, int hidden);
extern int dacq_hide_matching(char *name, int hidden);
extern char **dacq_get_lined_variables(int max_line);
extern int dacq_load_settings(char *filename);
extern int dacq_setup_variables(char *file);
/* scan.c */
extern int dacq_set_scan(int *chans, int nchan, double sfreq);
/* ch_convert.c */
extern int dacq_rack2scan(int rackCh);
extern int dacq_scan2rack(int scanCh);
extern int dacq_rack2scan_mcg(int rackCh);
extern int dacq_scan2rack_mcg(int scanCh);
extern int dacq_rack2log(int rackCh);
extern int dacq_log2rack(int logCh);
/* semaphores.c */
extern int dacq_check_sem(int id);
extern int dacq_get_sem(int id);
extern int dacq_release_sem(int id);
extern void dacq_setup_sems(void);
extern void dacq_remove_sems(void);
/* log.c */
extern void dacq_log(char *format, ...);
extern void dacq_log_buf(char *format, ...);
extern int dacq_log_buf_pending(void);
extern void dacq_log_set_name(char *new_name);
extern void dacq_log_set_time(int new_time);
extern void dacq_perror(char *s);
/* shmem.c */
extern dacqShmBlock dacq_get_shmem_block(int shmem_buf);
extern dacqShmBlock dacq_get_shmem(void);
extern int dacq_release_shmem(void);
extern void dacq_setup_shmem(void);
extern void dacq_remove_shmem(void);
/* nm_libversion.c */
extern const char *lib_dacq_vc_string(void);
extern int lib_dacq_version_major(void);
extern int lib_dacq_version_minor(void);
extern int lib_dacq_patch_level(void);
extern const char *lib_dacq_date_string(void);

#if defined(__cplusplus)
}
#endif
#endif /* of file */

