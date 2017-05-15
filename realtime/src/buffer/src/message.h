/*
 * Copyright (C) 2008, Robert Oostenveld & Christian Hesse
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#ifndef MESSAGE_H
#define MESSAGE_H

#include "platform_includes.h"

#ifdef __cplusplus
extern "C" {
#endif

/* FIXME these are obvious at the moment, but should be formally defined */
typedef char      CHAR_T;
typedef float  FLOAT32_T;
typedef double FLOAT64_T;

/* the following types should be according to "ISO C99: 7.18 Integer types" (see /usr/include/stdint.h on OSX and Linux) */

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

/* the same version number in the "other" endianness */
#define VERSION_OE  (UINT16_T) (((VERSION & 0x00FF) << 8) | ((VERSION & 0xFF00) >> 8))

/* these define the commands that can be used, which are split over the two available bytes */
#define PUT_HDR    (UINT16_T)0x0101 /* decimal 257 */
#define PUT_DAT    (UINT16_T)0x0102 /* decimal 258 */
#define PUT_EVT    (UINT16_T)0x0103 /* decimal 259 */
#define PUT_OK     (UINT16_T)0x0104 /* decimal 260 */
#define PUT_ERR    (UINT16_T)0x0105 /* decimal 261 */

#define GET_HDR    (UINT16_T)0x0201 /* decimal 513 */
#define GET_DAT    (UINT16_T)0x0202 /* decimal 514 */
#define GET_EVT    (UINT16_T)0x0203 /* decimal 515 */
#define GET_OK     (UINT16_T)0x0204 /* decimal 516 */
#define GET_ERR    (UINT16_T)0x0205 /* decimal 517 */

#define FLUSH_HDR  (UINT16_T)0x0301 /* decimal 769 */
#define FLUSH_DAT  (UINT16_T)0x0302 /* decimal 770 */
#define FLUSH_EVT  (UINT16_T)0x0303 /* decimal 771 */
#define FLUSH_OK   (UINT16_T)0x0304 /* decimal 772 */
#define FLUSH_ERR  (UINT16_T)0x0305 /* decimal 773 */

#define WAIT_DAT   (UINT16_T)0x0402 /* decimal 1026 */
#define WAIT_OK    (UINT16_T)0x0404 /* decimal 1027 */
#define WAIT_ERR   (UINT16_T)0x0405 /* decimal 1028 */

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
/*
  this should never be used to put data into the buffer, but is handy for handling conversions of other data types
*/
#define DATATYPE_UNKNOWN (UINT32_T)0xFFFFFFFF

/* these are used in the specification of the event selection criteria */
#define EVENTSEL_TYPE   1
#define EVENTSEL_VALUE  2
#define EVENTSEL_SAMPLE 3     /* for an exact match */
#define EVENTSEL_MINSAMPLE 4
#define EVENTSEL_MAXSAMPLE 5

/*
  if event->def->sample == EVENT_AUTO_SAMPLE, automatically insert current sample index instead
*/
#define EVENT_AUTO_SAMPLE   -1

/*
  the following enumeration is for specifying types of chunks that may be present in the "buf" part of the Fieldtrip header.
*/
enum {
    /** FT_CHUNK_UNSPECIFIED refers to a binary blob of known length, but unknown contents.
        Clients encountering this can try to use auto-detection, or just ignore this chunk.
        Unknown chunk types should be treated in the same manner.           */
    FT_CHUNK_UNSPECIFIED = 0,

    /** FT_CHUNK_CHANNEL_NAMES contains the channel names in ASCII format. Each channel is
        represented as a 0-terminated string (includes the case of just a 0 for an empty string).
        Example: chunk_data = "Left\0Right\0" for stereo sound signals.     */
    FT_CHUNK_CHANNEL_NAMES = 1,

    /** FT_CHUNK_CHANNEL_FLAGS contains a 0-terminated string describing the type of flags,
        and after that N (=#channels) bytes describing each channel. This is useful for
        specifying that a channel can have a discrete number of different types, e.g.
        chunk_data = "meg_ad_eog\0\1\1\1\1\3\3\2\2" should be used for a system with 8 channels,
        the first four of which are for MEG, then 2 channels EOG, then 2 channels A/D.  */
    FT_CHUNK_CHANNEL_FLAGS = 2,

    /** FT_CHUNK_RESOLUTIONS contains N double precision values mapping from A/D values to physical
        quantities such as micro-Volts in EEG.      */
    FT_CHUNK_RESOLUTIONS = 3,

    /** FT_CHUNK_ASCII_KEYVAL contains an arbitrary number of key/value pairs, each of
        which is given as a 0-terminated string. An empty key (=double 0) indicates the
        end of the list. Example: "amplifier_gain\0high\0noise_reduction\0active\0\0".  */
    FT_CHUNK_ASCII_KEYVAL = 4,

    /** FT_CHUNK_NIFTI1 contains a NIFTI-1 header (348 bytes long) */
    FT_CHUNK_NIFTI1 = 5,

    /** FT_CHUNK_SIEMENS_AP contains Siemens Protocol data in ASCII format (string) */
    FT_CHUNK_SIEMENS_AP = 6,

    /** FT_CHUNK_CTF_RES4 contains a .res4 file as written by the CTF MEG acquisition software (binary) */
    FT_CHUNK_CTF_RES4 = 7,

    /** FT_CHUNK_NEUROMAG_HEADER contains a .fif file as written by the Neuromag MEG acquisition software (binary) */
    FT_CHUNK_NEUROMAG_HEADER = 8,

    /** FT_CHUNK_NEUROMAG_ISOTRAK contains a .fif file as written by the Neuromag MEG acquisition software (binary) */
    FT_CHUNK_NEUROMAG_ISOTRAK = 9,

    /** FT_CHUNK_NEUROMAG_HPIRESULT contains a .fif file as written by the Neuromag MEG acquisition software (binary) */
    FT_CHUNK_NEUROMAG_HPIRESULT = 10
};

#pragma pack(push,1)

/* a packet that is sent over the network (or to disk) should contain the following */
typedef struct {
    UINT16_T version;   /* see VERSION */
    UINT16_T command;   /* see PUT_xxx, GET_xxx and FLUSH_xxx */
    UINT32_T bufsize;   /* size of the buffer in bytes */
} messagedef_t;

/* the header definition is fixed, except for the channel labels */
typedef struct {
    UINT32_T  nchans;
    UINT32_T  nsamples;
    UINT32_T  nevents;
    FLOAT32_T fsample;
    UINT32_T  data_type;
    UINT32_T  bufsize;     /* size of the buffer in bytes */
} headerdef_t;

/* the data definition is fixed */
typedef struct {
    UINT32_T nchans;
    UINT32_T nsamples;
    UINT32_T data_type;
    UINT32_T bufsize;     /* size of the buffer in bytes */
} datadef_t;

/* the event definition is fixed */
typedef struct {
    UINT32_T type_type;   /* usual would be DATATYPE_CHAR */
    UINT32_T type_numel;  /* length of the type string */
    UINT32_T value_type;
    UINT32_T value_numel;
    INT32_T  sample;
    INT32_T  offset;
    INT32_T  duration;
    UINT32_T bufsize;     /* size of the buffer in bytes */
} eventdef_t;

typedef struct {
    messagedef_t *def;
    void         *buf;
} message_t;

typedef struct {
    headerdef_t *def;
    void        *buf;     /* this can contain additional chunks, e.g. with the channel names */
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
    UINT32_T begsample; /* indexing starts with 0, should be >=0 */
    UINT32_T endsample; /* indexing starts with 0, should be <header.nsamples */
} datasel_t;

typedef struct {
    UINT32_T begevent;
    UINT32_T endevent;
} eventsel_t;

typedef struct {
    UINT32_T nsamples;
    UINT32_T nevents;
} samples_events_t;

typedef struct {
    samples_events_t threshold;
    UINT32_T milliseconds;
} waitdef_t;

typedef struct {
    UINT32_T type;      /* One of FT_CHUNK_** (see above) */
    UINT32_T size;      /* Size of chunk.data, total size is given by adding sizeof(ft_chunkdef_t)=8 */
} ft_chunkdef_t;

typedef struct {
    ft_chunkdef_t def;  /* See above. Note that this is not a pointer! */
    char data[1];       /* Data contained in this chunk */
} ft_chunk_t;

#pragma pack(pop)

#ifdef __cplusplus
}
#endif

#endif
