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
#include "rdadefs.h"

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

/** 'select' cannot handle more than 64 elements on Windows, but this
	should really be enough for all practical purposes. Depending on
	the sampling rate and number of channels, you would probably hit
	other performance boundaries first. Since the server socket itself
	also needs listening to (taking 1 away from the available 64), 
	NEVER set the following number to more than	63!!!
*/
#define RDA_MAX_NUM_CLIENTS  32

/** Number of blocks any client can lag behind before being disconnected */
#define RDA_MAX_LAG 5
 
/** RDA server control structure for starting, inspecting, and stopping a server */
typedef struct {
	pthread_t thread;			/**< Thread handle */
	pthread_mutex_t mutex;		/**< Mutex for protecting num_clients (actually not really necessary) */
	SOCKET server_socket;		/**< The server socket that clients connect to */
	int ft_buffer;				/**< Connection to FieldTrip buffer (socket or 0 for dmarequests) */
	volatile int num_clients;	/**< Current number of clients */
	volatile int should_exit;	/**< Flag to notify the server thread that it should stop */
	volatile int is_running;	/**< Flag that indicates whether the thread is still running */
	int blocksize;				/**< Block size for streaming out samples, 0 => adapt to incoming data */
	int use16bit;				/**< Flag that indicates whether 16 bit data should be streamed */
	int verbosity;				/**< Option that determines how much status information is printed during operation */
} rda_server_ctrl_t;

/** Internally used data structure to keep a linked list of 
	data packets that need to be sent out */
typedef struct rda_buffer_item {
	void *data;						/**< Points to complete RDA packet */
	size_t size;					/**< Size of the packet (=allocated memory block) */
	int blockNumber;				/**< Number of this data block (or -1 for start packet) */
	unsigned int refCount;			/**< Reference count (multiple clients get the same data) */
	struct rda_buffer_item *next;	/**< Next item in list or NULL */
} rda_buffer_item_t;

/** Internally used data structure to describe a client and its pending jobs */
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
	@param blocksize    Block size for streaming out samples (0=send out variable blocks depending on incoming data)
	@param errval		Optional pointer to an integer error value. Will contain either
						FT_NO_ERROR, FT_ERR_SOCKET, FT_OUT_OF_MEM or FT_THREADING on exit.
	@return 	Pointer to RDA server control structure, or NULL if an error occured
*/
rda_server_ctrl_t *rda_start_server(int ft_buffer, int use16bit, int port, int blocksize, int *errval);

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
