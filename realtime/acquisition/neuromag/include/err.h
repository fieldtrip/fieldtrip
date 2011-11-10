/* Automatically generated file, chages will not survive!*/
/* Copyright (c) 2007 Elekta Neuromag Oy, all rights reserved.*/
/* Generated  Fri Jun 15 15:36:16 EETDST 2007  by  skesti @ lupu.neuromag.fi  */

#ifndef err_h_included
#define err_h_included
/*
 * Header for library err.
 * 
 * Use this header for console programs. There is another library
 * called 'xerr.h' that contains the X-windows related calls that
 * are dependent on Motif headers.
 *
 * Copyright 2001 Neuromag Oy
 *
 * All rights reserved.
 * No part of this program may be photocopied, reproduced,
 * or translated to another program language without the
 * prior written consent of the author.
 *
 * $Header: err.hh,v 1.5 2004/09/30 08:34:16 mjk Exp $
 * $Id$
 */

#include <stdio.h>
#include <stdarg.h>

/**
 * Return values used by the library.   
 * \note failure code should be negative, to separate it from
 * valid stack depth return.
 */

#define ERR_RETURN_OK 0
#define ERR_RETURN_FAIL -1

/**
 * Limits.
 */

#define ERR_STACK_LEVELS_MAX  10
#define ERR_STRING_LEN_MAX    1024
#define ERR_MESSAGE_LEN_MAX   (1024*5)
#define ERR_MEMORY_SIZE       (ERR_STACK_LEVELS_MAX * ERR_STRING_LEN_MAX + ERR_MESSAGE_LEN_MAX + 36)

/**
 * Modes in err_get
 */

#define ERR_KEEP 1
#define ERR_FLUSH 0
#if defined(__cplusplus)
extern "C" {
#endif
/* errors.c */
extern void err_clear(void);
extern int err_num(void);
extern int err_push(const char *estr);
extern int err_set(const char *estr);
extern int err_pushf(const char *format, ...);
extern int err_setf(const char *format, ...);
extern int err_push_system(void);
extern int err_set_system(void);
extern char *err_pop(void);
extern char *err_peek(int level);
extern char *err_get(const char *epre, const char *esep, const char *epost, int keep);
extern int err_get_discarded(void);
extern void err_clear_error(void);
extern void err_set_error(const char *estr);
extern void err_set_sys_error(const char *estr);
extern void err_printf_set_error(const char *format, ...);
extern char *err_get_error(void);
extern void err_print_error(void);
/* err_single.c */
extern int err_init_single(void);
/* err_thread.c */
extern int err_init_thread(char *memory, int memsize);
/* nm_libversion.c */
extern const char *lib_err_vc_string(void);
extern int lib_err_version_major(void);
extern int lib_err_version_minor(void);
extern int lib_err_patch_level(void);
extern const char *lib_err_date_string(void);

#if defined(__cplusplus)
}
#endif
#endif /* of file */

