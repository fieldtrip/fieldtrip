/* 
 * Copyright (C) 2008, Robert Oostenveld & Christian Hesse
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: message.h,v $
 * Revision 1.11  2008/06/19 20:33:01  roboos
 * switched eventsel to signed integers, otherwise too difficult behaviour if they become negative
 *
 * Revision 1.10  2008/04/15 14:08:57  thohar
 * changed types of datasel_t.begsample and endsample to INT32_T to make it possible to set endsample=-1
 *
 * Revision 1.9  2008/04/14 14:11:54  thohar
 * added extern "C" statement for c++ builds
 *
 * Revision 1.8  2008/03/10 09:44:55  roboos
 * removed headersel, added property
 *
 * Revision 1.7  2008/02/26 21:43:26  roboos
 * renamed packet_t structure definition into messagedef_t, added message_t structure (contains def+buf)
 *
 * Revision 1.6  2008/02/20 13:36:03  roboos
 * changed comments to ansi style, needed for matlab
 *
 * Revision 1.5  2008/02/19 10:24:26  roboos
 * added copyright statement
 *
 * Revision 1.4  2008/02/18 17:04:08  roboos
 * lots of small changes, debugging for brainamp
 *
 * Revision 1.3  2008/02/18 12:57:44  roboos
 * fixed bug in typedefs for 32/64 bit linux versions
 *
 * Revision 1.2  2008/02/18 12:13:46  roboos
 * moved executable from buffer to demo
 * fixed bugs in sinewave and socket for events
 * stripped down the eventdef_t fields
 * many small changes
 *
 * Revision 1.1  2008/02/18 10:05:26  roboos
 * restructured the directory layout, copied all code into src, added directory for external code
 *
 * Revision 1.10  2008/02/14 16:57:46  roboos
 * removed old dcumentation
 * stripped down the content of eventdef_t
 * added propertydef_t
 * implemented eventsel_t with begin and endevent, c.f. datasel_t
 *
 * Revision 1.9  2008/02/11 21:40:53  roboos
 * fixed uint32 and 64, otherwise the wordsize check in buffer.c would fail on my powerbook
 *
 * Revision 1.8  2008/02/11 09:42:47  roboos
 * removed ; from end of line in WORDSIZE defines
 *
 * Revision 1.7  2008/02/11 09:32:22  roboos
 * removed all stuff that is not needed to get first implementation up and running
 * renamed structures (removed v1_)
 * added wrapper structures for headerdef+headerbuf instead of having a void * in headerdef
 * idem datadef and eventdef
 * be explicit about the uint16 type of version, command and datatypes
 *
 * Revision 1.6  2008/02/09 13:38:55  chrhes
 * corrected syntax errors in the definition of custom integer types that  had
 * caused UINT32_T and INT32_T to be defined as 64 bit integers
 *
 * Revision 1.5  2008/02/06 10:44:57  chrhes
 * added revision comment log
 *
 */

/* prevent double include */
#ifndef MESSAGE_H
#define MESSAGE_H

#include <stdint.h>
#ifdef __cplusplus
extern "C" {
	
#endif

/* FIXME these are obvious at the moment, but should be formally defined */
typedef char      CHAR_T;
typedef float  FLOAT32_T;
typedef double FLOAT64_T;

/* the following types should be according to "ISO C99: 7.18 Integer types" (see /usr/include/stdint.h on OSX and Linux) */
/* FIXME different endianness between client/server is not supported at the moment */

#ifndef INT8_T
typedef   int8_t   INT8_T;
#endif

#ifndef INT16_T
typedef  int16_t  INT16_T;
#endif

#ifndef INT32_T
typedef  int32_t  INT32_T;
#endif

#ifndef INT64_T
typedef  int64_t  INT64_T;
#endif

#ifndef UINT8_T
typedef  uint8_t  UINT8_T;
#endif

#ifndef UINT16_T
typedef uint16_t UINT16_T;
#endif

#ifndef UINT32_T
typedef uint32_t UINT32_T;
#endif

#ifndef UINT64_T
typedef uint64_t UINT64_T;
#endif

/* these can be used for indexing the buffer pointer as array */
#define WORDSIZE_CHAR    sizeof(CHAR_T   )
#define WORDSIZE_UINT8   sizeof(UINT8_T  ) 
#define WORDSIZE_UINT16  sizeof(UINT16_T )
#define WORDSIZE_UINT32  sizeof(UINT32_T )
#define WORDSIZE_UINT64  sizeof(UINT64_T )
#define WORDSIZE_INT8    sizeof(INT8_T   )
#define WORDSIZE_INT16   sizeof(INT16_T  )
#define WORDSIZE_INT32   sizeof(INT32_T  )
#define WORDSIZE_INT64   sizeof(INT64_T  )
#define WORDSIZE_FLOAT32 sizeof(FLOAT32_T)
#define WORDSIZE_FLOAT64 sizeof(FLOAT64_T)

/* define the version of the message packet */
#define VERSION    (UINT16_T)0x0001

/* these define the commands that can be used, which are split over the two available bytes */
#define PUT_HDR    (UINT16_T)0x0101
#define PUT_DAT    (UINT16_T)0x0102
#define PUT_EVT    (UINT16_T)0x0103
#define PUT_OK     (UINT16_T)0x0104
#define PUT_ERR    (UINT16_T)0x0105
#define PUT_PRP    (UINT16_T)0x0106

#define GET_HDR    (UINT16_T)0x0201 
#define GET_DAT    (UINT16_T)0x0202
#define GET_EVT    (UINT16_T)0x0203
#define GET_OK     (UINT16_T)0x0204
#define GET_ERR    (UINT16_T)0x0205
#define GET_PRP    (UINT16_T)0x0206

#define FLUSH_HDR  (UINT16_T)0x0301 
#define FLUSH_DAT  (UINT16_T)0x0302
#define FLUSH_EVT  (UINT16_T)0x0303
#define FLUSH_OK   (UINT16_T)0x0304
#define FLUSH_ERR  (UINT16_T)0x0305
#define FLUSH_PRP  (UINT16_T)0x0306

/* these are used in the data_t and event_t structure */
#define DATATYPE_CHAR    (UINT32_T)0
#define DATATYPE_UINT8   (UINT32_T)1
#define DATATYPE_UINT16  (UINT32_T)2
#define DATATYPE_UINT32  (UINT32_T)3
#define DATATYPE_UINT64  (UINT32_T)4
#define DATATYPE_INT8    (UINT32_T)5
#define DATATYPE_INT16   (UINT32_T)6
#define DATATYPE_INT32   (UINT32_T)7
#define DATATYPE_INT64   (UINT32_T)8
#define DATATYPE_FLOAT32 (UINT32_T)9
#define DATATYPE_FLOAT64 (UINT32_T)10

/* these are used in the specification of the event selection criteria */
#define EVENTSEL_TYPE   1
#define EVENTSEL_VALUE  2
#define EVENTSEL_SAMPLE 3     /* for an exact match */
#define EVENTSEL_MINSAMPLE 4
#define EVENTSEL_MAXSAMPLE 5

/* a packet that is sent over the network (or to disk) should contain the following */
typedef struct {
	UINT16_T version;   /* see VERSION */
	UINT16_T command;   /* see PUT_xxx, GET_xxx and FLUSH_xxx */
	UINT32_T bufsize;   /* size of the buffer in bytes */
} messagedef_t;

/* the header definition is fixed, except for the channel labels */
typedef struct {
	UINT32_T nchans;
	UINT32_T nsamples;
	UINT32_T nevents;
	FLOAT32_T fsample;
	UINT32_T data_type;
	UINT32_T bufsize;   /* size of the buffer in bytes */
} headerdef_t;

/* the data definition is fixed */
typedef struct {
	UINT32_T nchans;
	UINT32_T nsamples;
	UINT32_T data_type;
	UINT32_T bufsize;   /* size of the buffer in bytes */
} datadef_t;

/* the event definition is fixed */
typedef struct {
	UINT32_T type_type;        /* usual would be DATATYPE_CHAR */
	UINT32_T type_numel;       /* length of the type string */
	UINT32_T value_type;
	UINT32_T value_numel; 
	INT32_T sample; 
	INT32_T offset;
	INT32_T duration;
	UINT32_T bufsize;          /* size of the buffer in bytes */
} eventdef_t;

typedef struct {
	UINT32_T type_type;        /* usual would be DATATYPE_CHAR */
	UINT32_T type_numel;       /* length of the type string */
	UINT32_T value_type;
	UINT32_T value_numel; 
	UINT32_T bufsize;          /* size of the buffer in bytes */
} propertydef_t;

typedef struct {
	messagedef_t *def;
	void         *buf;
} message_t;

typedef struct {
	headerdef_t *def;
	void        *buf;			/* FIXME this should contain the channel names */
} header_t;

typedef struct {
	datadef_t  *def;
	void       *buf;
} data_t;

typedef struct {
	eventdef_t *def;
	void       *buf;
} event_t;

typedef struct {
	propertydef_t *def;
	void          *buf;
} property_t;

typedef struct {
	INT32_T begsample; /* indexing starts with 0, should be >=0 */
	INT32_T endsample; /* indexing starts with 0, should be <header.nsamples */
} datasel_t;

typedef struct {
	INT32_T begevent;
	INT32_T endevent;
} eventsel_t;

#ifdef __cplusplus
}
#endif

#endif
