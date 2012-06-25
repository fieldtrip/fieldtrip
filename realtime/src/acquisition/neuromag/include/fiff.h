/* Automatically generated file, chages will not survive!*/
/* Copyright (c) 2010 Elekta Neuromag Oy, all rights reserved.*/
/* Generated  Mon Apr 19 11:22:22 EEST 2010  by  mjk @ mjk.neuromag.fi  */

#ifndef fiff_h_included
#define fiff_h_included
/*
 *
 * Copyright 1990-2003 Elekta Neuromag Oy
 *
 * No part of this program may be copied, photocopied, reproduced,
 * or translated to another program language without the
 * prior written consent of the author.
 *
 * $Header: /sw/tmp/build/libs-fiff-2.1.0/src/RCS/fiff.hh,v 1.7 2006/04/26 11:35:08 mjk Exp $
 *
 * This is a header file for library fiff, providing access to 
 * low level functions that enable reading and writing of FIFF files.
 *
 */

#include <stdio.h>
#include "fiff_file.h"
#include "fiff_types.h"


#if defined(__cplusplus)
extern "C" {
#endif
/* fiff_compat.c */
extern void fiff_convert_shorts(short *data, size_t ndata);
extern void fiff_convert_ints(int *data, size_t ndata);
extern void fiff_convert_floats(float *data, size_t ndata);
extern void fiff_convert_doubles(double *data, size_t ndata);
extern void fiff_convert_shorts_to_floats(short *data, size_t ndata);
extern void fiff_convert_ints_to_floats(int *data, size_t ndata);
extern void fiff_convert_tag_data(fiffTag tag, int from);
extern void fiff_convert_tag_info(fiffTag tag);
extern int fiff_get_time(long *secs, long *usecs);
extern int fiff_protect_file(const char *path);
extern int fiff_deprotect_file(const char *path);
extern int fiff_truncate_file(int fildes, off_t length);
/* fiff_io.c */
extern int fiff_safe_fseek(FILE *fp, long offset, int whence);
extern int fiff_safe_fread(void *buf, size_t size, FILE *fp);
extern int fiff_safe_fwrite(void *buf, size_t size, FILE *fp);
extern int fiff_read_tag_header(FILE *in, long pos, fiffTag tag);
extern int fiff_read_tag_data(FILE *in, fiffTag tag);
extern int fiff_advance_after_header(FILE *in, fiffTag tag);
extern int fiff_advance_after_data(FILE *in, fiffTag tag);
extern int fiff_read_tag_info(FILE *in, fiffTag tag);
extern int fiff_read_this_tag(FILE *in, long lpos, fiffTag tag);
extern int fiff_read_tag(FILE *in, fiffTag tag);
extern int fiff_load_this_buffer(FILE *in, long lpos, fiffTag tag);
extern int fiff_load_buffer(FILE *in, fiffTag tag);
extern int fiff_write_tag_header(FILE *file, long pos, fiffTag tag);
extern int fiff_write_tag_data(FILE *file, fiffTag tag);
extern int fiff_write_tag_info(FILE *file, fiffTag tag);
extern int fiff_write_this_tag(FILE *file, long lpos, fiffTag tag);
extern int fiff_write_tag(FILE *file, fiffTag tag);
extern int fiff_start_block(FILE *file, int kind);
extern int fiff_end_block(FILE *file, int kind);
extern int fiff_start_file(FILE *file);
extern int fiff_end_file(FILE *file);
extern int fiff_write_floats(FILE *file, float *data, size_t ndata);
extern int fiff_write_doubles(FILE *file, double *data, size_t ndata);
extern int fiff_write_ints(FILE *file, int *data, size_t ndata);
/* fiff_id.c */
extern int fiff_new_file_id(fiffId id);
extern int fiff_get_date(long *date);
extern int fiff_id_match(fiffId id1, fiffId id2);
/* fiff_explain.c */
extern void fiff_explain(int kind);
extern char *fiff_get_tag_explanation(int kind);
extern void fiff_explain_block(int kind);
extern char *fiff_get_block_explanation(int kind);
extern char *fiff_explain_unit(int unit, int mul);
/* fiff_dir.c */
extern int fiff_make_dir(fiffFile file, int ptrpos, int dirpos, int make_dir);
extern int fiff_read_dir(fiffFile file, int ptrpos, int dirpos, int make_dir);
extern int fiff_how_many_entries(fiffDirEntry dir);
extern int fiff_put_dir(FILE *fd, fiffDirEntry dir);
/* fiff_open.c */
extern void fiff_close(fiffFile file);
extern fiffFile fiff_open_quick(const char *name);
extern fiffFile fiff_open(const char *name);
extern fiffFile fiff_open_update(const char *name);
extern fiffFile fiff_open_rescue(const char *name);
extern fiffFile fiff_open_fix(const char *name);
/* fiff_julian.c */
extern int fiff_julday(int id, int mm, int iyyy, long *julian);
extern int fiff_caldate(long julian, int *id, int *mm, int *iyyy);
/* fiff_insert.c */
extern int fiff_insert_after(fiffFile dest, int where, fiffTag tags, int ntag);
/* fiff_trans.c */
extern fiffCoordTrans fiff_make_transform(int from, int to, float rot[3][3], float move[3]);
extern void fiff_coord_trans(float r[3], fiffCoordTrans t, int do_move);
extern void fiff_coord_trans_inv(float r[3], fiffCoordTrans t, int do_move);
extern fiffCoordTrans fiff_make_transform_card(int from, int to, float *rL, float *rN, float *rR);
extern fiffCoordTrans fiff_invert_transform(fiffCoordTrans t);
extern fiffCoordTrans fiff_dup_transform(fiffCoordTrans t);
extern fiffCoordTrans fiff_combine_transforms(int from, int to, fiffCoordTrans t1, fiffCoordTrans t2);
/* fiff_writes.c */
extern int fiff_write_string_tag(FILE *out, int kind, const char *text);
extern int fiff_write_int_tag(FILE *out, int kind, fiff_int_t val);
extern int fiff_write_float_tag(FILE *out, int kind, double val);
extern int fiff_write_coord_trans(FILE *out, fiffCoordTrans t);
extern fiffId fiff_write_id(FILE *out, int kind);
extern int fiff_write_this_id(FILE *out, int kind, fiffId id);
/* fiff_def_dir.c */
extern char *fiff_def_name_dir(const char *name);
/* fiff_dir_tree.c */
extern void fiff_dir_tree_free(fiffDirNode node);
extern int fiff_dir_tree_create(fiffFile file, fiff_data_t **dat);
extern void fiff_dir_tree_print(fiffDirNode tree);
extern int fiff_dir_tree_count(fiffDirNode tree);
extern fiffDirNode *fiff_dir_tree_find(fiffDirNode tree, int kind);
extern fiffDirNode fiff_dir_tree_find_id(fiffDirNode tree, fiffId id);
extern fiffTag fiff_dir_tree_get_tag(fiffFile file, fiffDirNode node, int kind);
extern int fiff_copy_subtree(FILE *to, fiffFile from, fiffDirNode node);
/* fiff_matrix.c */
extern int *fiff_get_matrix_dims(fiffTag tag);
extern float **fiff_get_float_matrix(fiffTag tag);
extern double **fiff_get_double_matrix(fiffTag tag);
extern int **fiff_get_int_matrix(fiffTag tag);
extern int fiff_write_float_matrix(FILE *out, int kind, float **data, int rows, int cols);
extern int fiff_write_double_matrix(FILE *out, int kind, double **data, int rows, int cols);
extern int fiff_write_int_matrix(FILE *out, int kind, int **data, int rows, int cols);
/* fiff_pack.c */
extern void fiff_pack_data(const float *data, short *packed_data, float *offset, float *scale, int nsamp);
extern void fiff_unpack_data(double offset, double scale, const short *packed, int nsamp, float *orig);
extern void fiff_unpack_HP3852A(const short *packed, int nrd, float *unpacked);
/* fiff_type_spec.c */
extern fiff_int_t fiff_type_fundamental(fiff_int_t type);
extern fiff_int_t fiff_type_base(fiff_int_t type);
extern fiff_int_t fiff_type_matrix_coding(fiff_int_t type);
/* fiff_sparse.c */
extern void fiff_free_sparse_matrix(fiffSparseMatrix m);
extern fiff_int_t *fiff_get_matrix_sparse_dims(fiffTag tag);
extern fiffSparseMatrix fiff_get_float_sparse_matrix(fiffTag tag);
extern int fiff_write_float_sparse_matrix(FILE *out, int kind, fiffSparseMatrix mat);
/* fiff_subsystem.c */
extern void fiff_free_hpi_subsys(fiffHpiSubsys hpi);
extern fiffHpiSubsys fiff_get_hpi_subsys(fiffFile file, fiffDirNode node);
/* nm_libversion.c */
extern const char *lib_fiff_vc_string(void);
extern int lib_fiff_version_major(void);
extern int lib_fiff_version_minor(void);
extern int lib_fiff_patch_level(void);
extern const char *lib_fiff_date_string(void);

#if defined(__cplusplus)
}
#endif
#endif /* of file */

