/*
 * Copyright (C) 2010, Stefan Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#ifndef __SOCKETSERVER_H
#define __SOCKETSERVER_H

#include <pthread.h>
#include "buffer.h"

#define MERGE_THRESHOLD 4096    /* TODO: optimize this value? Maybe look at MTU size */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef WIN32
#define INVALID_SOCKET  -1
typedef int SOCKET;
#include <sys/un.h>
#else
typedef int socklen_t;
#endif

typedef int (*ft_request_callback_t)(const message_t *request, message_t **response, void *user_data);

/** The following structure is used for managing a server. The structure is
    allocated and filled in ft_start_buffer_server and then passed on to the
    actual server thread, as well as to all threads handling the client 
    connections. Thread functions are written such that they monitor the
    "keepRunning" member of this structure, and once this is set to 0, the
    threads stop and exit.
*/
typedef struct {
        SOCKET serverSocket;            /**< Socket the server listens on (TCP or UNIX domain) */
        int keepRunning;                /**< Flag=1: thread functions keep looping, 0: threads exit */
        int numClients;                 /**< Current number of clients connected to the server */
        int verbosity;                  /**< Determines how much information is being printed during operation */
        int isUnixDomain;               /**< 1: UNIX domain socket, 0: TCP socket */
        pthread_t threadID;             /**< POSIX thread identifier of the server thread, client threads are detached immediately */
        pthread_mutex_t lock;           /**< Mutex to protect the "numClients" member, commonly used by all threads */
        ft_request_callback_t callback; /**< Callback function to be called *instead* of dmarequest */
        void *user_data;                /**< Pointer to user-defined data structure, passed on to callback */
} ft_buffer_server_t;

/** Small helper structure that is passed to client threads. Get's allocated
        using malloc() in the server thread, and disposed using free() in the client
        thread.
*/
typedef struct {
        ft_buffer_server_t *server;     /**< Pointer to the common control structure */
        SOCKET clientSocket;            /**< The newly created socket (from "accept") */
        int mergePackets;               /**< 1: merge packets if total size below threshold, 0: never merge (=>local host) */
} ft_buffer_socket_t;

/** Creates a server socket, binds it to the specified UNIX domain name or port,
        starts listening on this socket, and spawns a background thread to serve
        requests on this socket.
        TCP sockets are created for positive port numbers (name is ignored), 
        UNIX domain sockets are created for port=0 (name needs to be a UNIX pathname).

        If callback!=NULL, that function is called like
                callback(request, &response, user_data)
        for every request coming in over the socket, *instead* of the normal dmarequest.
        Usually the user will call the latter function internally. Make sure your
        callback is re-entrant!!!
        
        Returns allocated control structure, or NULL in case of errors.
*/      
ft_buffer_server_t *ft_start_buffer_server(int port, const char *name, ft_request_callback_t callback, void *user_data);

/** Stops background thread(s), closes the sockets, and disposes the control structure S.
        S cannot be used anymore after this call. 
*/
void ft_stop_buffer_server(ft_buffer_server_t *S);

#ifdef __cplusplus
}
#endif

#endif /* __SOCKETSERVER_H */
