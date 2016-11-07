/*
 * Collection of routines for saving FieldTrip buffer data to disk.
 *
 * (C) 2010 S. Klanke
 */

#ifndef __ft_storage_h
#define __ft_storage_h

#include <stdio.h>
#include <buffer.h>

#ifdef __cplusplus
extern "C" {
#endif

/* FIXME: move this somewhere else */
#define FT_OUT_OF_MEMORY 1000
#define FT_FILE_ERROR    1001

typedef struct {
	int created;
	FILE *fHeader;
	FILE *fHeaderTxt;
	FILE *fSamples;
	FILE *fEvents;
	FILE *fTime;
	char *dirName;
	int dirLen;
	int curSampleFile;
	UINT32_T sampleSize; /* wordsize(data_type) * numChannels */
	UINT32_T numChannels;
	UINT32_T numSamples;
	UINT32_T numEvents;
} ft_storage_t;

typedef struct {
	double time;
	int numSamples;
	int numEvents;
} ft_timing_element_t;

ft_storage_t *ft_storage_create(const char *directory, const headerdef_t *hdef, const void *chunks, int *errCode);
void ft_storage_close       (ft_storage_t *S);
int  ft_storage_add_samples (ft_storage_t *S, int numSamples, const void *data);
int  ft_storage_add_events  (ft_storage_t *S, int size, const void *events);
int  ft_storage_add_event   (ft_storage_t *S, const eventdef_t *event, const void *type, const void *value);
int  ft_storage_add_timing  (ft_storage_t *S, const ft_timing_element_t *te);

/*
ft_storage_t *ft_storage_open(const char *directory);
headerdef_t *ft_storage_read_header(const ft_storage_t *S);
int ft_storage_read_data(ft_storage_t *S, int begsample, int endsample, void *data);
void *ft_storage_read_events(ft_storage_t *S);
*/

#ifdef __cplusplus
}
#endif

#endif
