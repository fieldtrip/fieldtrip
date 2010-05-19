/*
 * Copyright (C) 2010 S. Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 *
 */

#ifndef __rdaserver_h
#define __rdaserver_h
 
#include <pthread.h>
#include "buffer.h" 

#ifdef __cplusplus
extern "C" {
#endif

/* On Windows, sockets are not described by a plain int, but by the type SOCKET, which has
	the size of a pointer. On WIN64, sizeof(SOCKET) != sizeof(int), although some people still
	argue that it is safe to cast between those. To be sure, the RDA server implementation 
	always uses SOCKET as the base type, and defines this as an 'int' on POSIX systems.
*/
#ifndef WIN32
typedef int SOCKET;
#define INVALID_SOCKET  -1
#endif

/** Error values as returned by rda_start_server */
#define FT_NO_ERROR         0
#define FT_ERR_OUT_OF_MEM   1
#define FT_ERR_SOCKET       2
#define FT_ERR_THREADING    3

/** Message types as sent to RDA clients */
#define RDA_START_MSG	1
#define RDA_INT_MSG		2
#define RDA_STOP_MSG	3
#define RDA_FLOAT_MSG	4
 
/** RDA server control structure for starting, inspecting, and stopping a server */
typedef struct {
	pthread_t thread;			/**< Thread handle */
	pthread_mutex_t mutex;		/**< Mutex for protecting num_clients (actually not really necessary) */
	SOCKET server_socket;		/**< The server socket that clients connect to */
	int ft_buffer;				/**< Connection to FieldTrip buffer (socket or 0 for dmarequests) */
	volatile int num_clients;	/**< Current number of clients */
	volatile int should_exit;	/**< Flag to notify the server thread that it should stop */
	volatile int is_running;	/**< Flag that indicates whether the thread is still running */
	int use16bit;				/**< Flag that indicates whether 16 bit data should be streamed */
	int verbosity;				/**< Option that determines how much status information is printed during operation */
} rda_server_ctrl_t;

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

/** Internally used data structure to keep a linked list of 
	data packets that need to be sent out */
typedef struct rda_buffer_item {
	void *data;						/**< Points to complete RDA packet */
	size_t size, sizeAlloc;			/**< Size of the packet and size of allocated memory block */
	unsigned int refCount;			/**< Reference count (multiple clients get the same data) */
	struct rda_buffer_item *next;	/**< Next item in list or NULL */
} rda_buffer_item_t;

/** Internally used data structure to describe a bunch of clients
	and their pending jobs */
typedef struct {
	SOCKET sock;				/**< Client socket */
	rda_buffer_item_t *item;	/**< Points to the current/next packet to be written */
	size_t written;				/**< Number of bytes that have been written (from item->data) */
} rda_client_job_t;

/** Helper function for converting any FieldTrip data type to single precision floats
	@param N			number of values to convert
	@param dest 		destination, must point to an array of at least N floats
	@param data_type	data type as described by FieldTrip DATATYPE_** constants
	@param src			source buffer
*/
void rda_aux_convert_to_float(UINT32_T N, void *dest, UINT32_T data_type, const void *src);

/** Starts an RDA server with a given FieldTrip connection (usually 0 for DMA), serving
    either single precision or 16 bit integer data.
	@param ft_buffer	FieldTrip connection (0 for DMA, or socket for TCP connection)
	@param use16bit		pass 0 to serve single precision data (with conversion as necessary)
						pass non-zero to serve 16 bit integers (only works if the FieldTrip buffer
						also contains 16 bit data)
	@param port			Port number to bind to, or 0 for default port (51244 for 16 bit, 51344 for single precision)
	@param errval		Optional pointer to an integer error value. Will contain either
						FT_NO_ERROR, FT_ERR_SOCKET, FT_OUT_OF_MEM or FT_THREADING on exit.
	@return 	Pointer to RDA server control structure, or NULL if an error occured
*/
rda_server_ctrl_t *rda_start_server(int ft_buffer, int use16bit, int port, int *errval);

/** Stops an RDA server by closing the associated connections, stopping the background thread, and
	deallocating its memory (including the control structure pointed to by the argument)
	@param	SC		Must point to an RDA server control structure as created by rda_start_server
	@return -1	if SC is NULL
			0	on success
*/
int rda_stop_server(rda_server_ctrl_t *SC);

/** Simple helper function for retrieving the current number of clients of an RDA server 
	@param	SC		Must point to an RDA server control structure as created by rda_start_server
	@return	The number of clients
*/	
int rda_get_num_clients(rda_server_ctrl_t *SC);

#ifdef __cplusplus
}
#endif

#endif /* __rdaserver_h */
