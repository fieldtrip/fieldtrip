/*
 * Copyright (C) 2010 S. Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 *
 */

#ifndef __rdadefs_h
#define __rdadefs_h

#include "message.h"		/* for the integer datytypes */
 
#ifdef __cplusplus
extern "C" {
#endif

/** Message types as sent to RDA clients */
#define RDA_START_MSG	1
#define RDA_INT_MSG		2
#define RDA_STOP_MSG	3
#define RDA_FLOAT_MSG	4

#pragma pack(push,1)

/** Structure of the first 24 bytes of all RDA messages */
typedef struct {
	UINT8_T guid[16];
	UINT32_T nSize; /* Size of the message block in bytes including this header */
	UINT32_T nType; /* 1:start 2:int16 3:stop 4:float */
} rda_msg_hdr_t;

/** Describes the structure of RDA data messages (fixed part only) */
typedef struct {
	rda_msg_hdr_t hdr;
	UINT32_T nBlock;	/* Block number, i.e. acquired blocks since acquisition started. */
	UINT32_T nPoints;	/* Number of data points (samples) in this block */
	UINT32_T nMarkers;  /* Number of markers in this block */
	/* after this, you get the data, and then an array of markers */
} rda_msg_data_t;

/** Describes the structure of an RDA start message (fixed part only) */
typedef struct {
	rda_msg_hdr_t hdr;
	UINT32_T nChannels;			/* Number of channels */
	double dSamplingInterval;	/* Sampling interval in microseconds */
	/* after this, you have double dResolutions[] and 
	   channels names come after this as 0-terminated strings */
} rda_msg_start_t;

/* TODO: we never send a STOP packet, because it's hard to tell when the buffer will
	stop receiving data. We should maybe try to detect if a new header is put into the
	buffer, and send a STOP and a START in that case.
*/	

/** Describes the structure of an RDA marker (fixed part only) */
typedef struct {
	UINT32_T nSize;		/* Size of this marker */
	UINT32_T nPosition; /* Relative position in the data block */
	UINT32_T nPoints;	/* Number of points of this marker */
	INT32_T nChannel;	/* Associated channel number (-1 = all channels) */
	/* char sTypeDesc[1];  Type description in ASCII delimited by '\0', variable length actually */
} rda_marker_t;	

#pragma pack(pop)

#ifdef __cplusplus
}
#endif

#endif /* __rdadefs_h */
