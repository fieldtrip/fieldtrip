/*
 * Copyright (C) 2010, Stefan Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fcntl.h>
#include <errno.h>
#include <socketserver.h>

/* this function deals with the incoming client request */
void *_buffer_socket_func(void *arg) {
	SOCKET sock;
	ft_buffer_server_t *SC;
	int mergePackets;
	messagedef_t reqdef;
	message_t request;
	message_t *response = NULL;
	int state=0; /* 0 = reading def, 1=reading buf, 2=writing def, 3=writing buf */
	int bytesDone, bytesTotal;
	char mergeBuffer[MERGE_THRESHOLD];
	char *curPtr;	/* points at buffer that needs to be filled or written out */
	int swap = 0;	
	int canRead, canWrite;
	UINT16_T reqCommand;
	UINT32_T respBufSize;
	struct timeval tv;
	fd_set readSet, writeSet;

	if (arg==NULL) return NULL;
		
	/* copy over necessary variables and free the given structure */
	SC 			 = ((ft_buffer_socket_t *) arg)->server;
	sock 		 = ((ft_buffer_socket_t *) arg)->clientSocket;
	mergePackets = ((ft_buffer_socket_t *) arg)->mergePackets;
	free(arg);
	
	if (SC->verbosity > 0) {
		printf("Started new client thread with packet merging = %i\n", mergePackets);
	}
	
    pthread_mutex_lock(&SC->lock);
    SC->numClients++;
    pthread_mutex_unlock(&SC->lock);
	
	request.def = &reqdef;
	request.buf = NULL;
	bytesDone = 0;
	bytesTotal = sizeof(messagedef_t);
	curPtr = (char *) request.def;
	
	tv.tv_sec = 0;
	tv.tv_usec = 10000;
	

	/* keep processing messages untill the connection is closed */
	while (SC->keepRunning) {
		int sel, res, n;

		FD_ZERO(&readSet);
		FD_ZERO(&writeSet);
		if (state < 2) {
			FD_SET(sock, &readSet);
		} else {
			FD_SET(sock, &writeSet);
		}
		sel = select(sock+1, &readSet, &writeSet, NULL, &tv);
		if (sel == 0) continue;
		if (sel < 0) {
			fprintf(stderr, "Error in 'select' operation - closing client connection.\n");
			break;
		}
		canRead = FD_ISSET(sock, &readSet);
		canWrite = FD_ISSET(sock, &writeSet);
		
		if (canRead) {
			n = recv(sock, curPtr + bytesDone, bytesTotal - bytesDone, 0);
			if (n<=0) {
				/* socket was closed */
				if (SC->verbosity>0) {
					printf("Remote side closed client connection\n");
				}
				break;
			}
			bytesDone+=n;
			if (bytesDone<bytesTotal) continue;
			
			if (state == 0) {
				/* we've read the request.def completely */
				if (reqdef.version==VERSION_OE) {
					swap = 1;
					ft_swap16(2, &reqdef.version); /* version + command */
					ft_swap32(1, &reqdef.bufsize);
					reqCommand = reqdef.command;		
				}
				if (reqdef.version!=VERSION) {
					fprintf(stderr,"Incorrect version requested - closing socket.\n");
					break;
				}
				if (reqdef.bufsize > 0) {
					request.buf = malloc(reqdef.bufsize);
					if (request.buf == NULL) {
						fprintf(stderr, "Out of memory\n");
						break;
					}
					curPtr = request.buf;
					bytesDone = 0;
					bytesTotal = reqdef.bufsize;
					state = 1;
					continue;
				} 
			} else {
				/* state must be 1, but we've read request.buf completely */	
				if (swap) ft_swap_buf_to_native(reqCommand, reqdef.bufsize, request.buf);
			}
			
			/* request read completely, deal with it */

			/* print_request(request.def);  */
			if (SC->callback != NULL) {
				res = SC->callback(&request, &response, SC->user_data);
				if (res != 0 || response == NULL || response->def == NULL) {
					fprintf(stderr, "buffer_socket_func: an unexpected error occurred in user-defined request handler\n");
					break;
				}
			} else {
				res = dmarequest(&request, &response);
				if (res != 0 || response == NULL || response->def == NULL) {
					fprintf(stderr, "buffer_socket_func: an unexpected error occurred in dmarequest\n");
					break;
				}
			}
		
			if (request.buf != NULL) {
				free(request.buf);
				request.buf = NULL;
			}
			
			respBufSize = response->def->bufsize;
			if (swap) ft_swap_from_native(reqCommand, response);
		
			/* merge response->def and response->buf if they are small, so we can send it in one go over TCP */
			if (mergePackets && respBufSize > 0 && respBufSize + sizeof(messagedef_t) <= MERGE_THRESHOLD) {
				memcpy(mergeBuffer, response->def, sizeof(messagedef_t));
				memcpy(mergeBuffer + sizeof(messagedef_t), response->buf, respBufSize);
				
				curPtr = mergeBuffer;
				bytesDone = 0;
				bytesTotal = respBufSize + sizeof(messagedef_t);
				state = 3;
			} else {
				curPtr = (char *) response->def;
				bytesDone = 0;
				bytesTotal = sizeof(messagedef_t);
				state = 2;
			}
			canWrite = 1;
		}
		
		if (state >= 2 && canWrite) {
			n = send(sock, curPtr + bytesDone, bytesTotal - bytesDone, 0);
			if (n<=0) {
				/* socket was closed */
				fprintf(stderr, "Cannot write to socket -- closing client connection.\n");
				break;
			}
			bytesDone+=n;
			if (bytesDone < bytesTotal) continue;
			if (state==2 && respBufSize > 0) {
				curPtr = (char *) response->buf;
				bytesDone = 0;
				bytesTotal = respBufSize;
				state = 3;
				continue;
			}
			/* done ! */
			if (response->buf) free(response->buf);
			free(response->def);
			free(response);
			response = NULL;
			state = 0;
			curPtr = (char *) request.def;
			bytesDone = 0;
			bytesTotal = sizeof(messagedef_t);
		}
	}
	
    pthread_mutex_lock(&SC->lock);
    SC->numClients--;
    pthread_mutex_unlock(&SC->lock);
	
	closesocket(sock);
	if (request.buf!=NULL) free(request.buf);
    if (response!=NULL) {
		if (response->buf!=NULL) free(response->buf);
		if (response->def!=NULL) free(response->def);
		free(response);
	}
	
	return NULL;
}



/***********************************************************************
 * this thread listens to incoming TCP/UNIX domain socket connections
 * if a connection is made by a client, it starts _buffer_socket_func
 ***********************************************************************/
void *_buffer_server_func(void *arg) {
	ft_buffer_server_t *SC = (ft_buffer_server_t *) arg;
	struct timeval tv;			/* for select timeout */
	int n;
	
	if (SC == NULL) {
		fprintf(stderr, "FieldTrip buffer server thread started with invalid argument\n");
		return NULL;
	}
	
	tv.tv_sec = 0;
	tv.tv_usec = 10000; /* 10 ms */
	
	while (SC->keepRunning) {
		SOCKET c;
		fd_set readSet;
		int sel, rc, merge;
		pthread_t tid;
		ft_buffer_socket_t *CC;
		
		/*
		 * If no pending connections are present on the queue, and the socket
		 * is not marked as non-blocking, accept() blocks the caller until a
		 * connection is present.  If the socket is marked non-blocking and
		 * no pending connections are present on the queue, accept() returns
		 * an error as described below.
		 */
		 
		FD_ZERO(&readSet);
		FD_SET(SC->serverSocket, &readSet);
		
		sel = select(SC->serverSocket + 1, &readSet, NULL, NULL, &tv);
		
		if (sel == 0) continue;
		
		/* The following code portion looks weird because of the preprocessor
			defines splitting an if/else clause, but it's just what we want:
			On Windows, there is no if/else clause, and the second (TCP) bit
			is always called. On POSIX systems, we branch depending on whether
			the server runs a local domain socket or TCP socket.
		*/
#ifndef WIN32		
		if (SC->isUnixDomain) {
			struct sockaddr_un sa;
			socklen_t size_sa = sizeof(sa);
			
			c = accept(SC->serverSocket, (struct sockaddr *)&sa, &size_sa);
			
			if (c == INVALID_SOCKET) {
				perror("buffer_server, accept");
				continue;
			}
			/* never merge packets for (local) UNIX sockets */
			merge = 0;
		} else 
#endif		
		{
			struct sockaddr_in sa;
			socklen_t size_sa = sizeof(sa);
			
			c = accept(SC->serverSocket, (struct sockaddr *)&sa, &size_sa);
			
			if (c == INVALID_SOCKET) {
				perror("buffer_server, accept");
				continue;
			}
			/* enable packet merging only if it's not localhost */
			merge = (sa.sin_addr.s_addr == htonl(INADDR_LOOPBACK)) ? 0 : 1;
		}
		
		CC = (ft_buffer_socket_t *) malloc(sizeof(ft_buffer_socket_t));
		if (CC==NULL) {
			fprintf(stderr, "Out of memory\n");
			continue;
		}
		
		CC->server = SC;
		CC->clientSocket = c;
		CC->mergePackets = merge;
		
		rc = pthread_create(&tid, NULL, _buffer_socket_func, CC);
		if (rc) {
			fprintf(stderr, "tcpserver: return code from pthread_create() is %d\n", rc);
			closesocket(c);
			free(CC);
		}
	}
	while ((n=SC->numClients)>0) {
		printf("Waiting for %i remaining client threads to stop...\n", n);
		usleep(10000);
	}
	pthread_mutex_lock(&SC->lock);
	SC->numClients--;
	pthread_mutex_unlock(&SC->lock);
	return NULL;
}


ft_buffer_server_t *ft_start_buffer_server(int port, const char *name, ft_request_callback_t callback, void *user_data) {
	ft_buffer_server_t *SC;
	int optval;
	SOCKET s = INVALID_SOCKET;
	
	SC = (ft_buffer_server_t *) malloc(sizeof(ft_buffer_server_t));
	if (SC==NULL) return NULL;
	
	SC->callback = callback;
	SC->user_data = user_data;

#ifdef WIN32
	{
		WSADATA wsa;
		if(WSAStartup(MAKEWORD(1, 1), &wsa))
		{
			fprintf(stderr, "ft_start_buffer_server: cannot start WIN32 sockets.\n");
			goto cleanup;
		}
	}
#endif

	/* setup socket */
	if (port == 0) {
#ifdef WIN32
		fprintf(stderr, "ft_start_buffer_server: invalid port number given.\n");
		goto cleanup;
#else
		struct sockaddr_un sa;
		/* UNIX domain socket */
		s = socket(AF_UNIX, SOCK_STREAM, 0);
		if (s == INVALID_SOCKET) {
			perror("ft_start_buffer_server, socket");
			goto cleanup;
		}
		sa.sun_family = AF_UNIX;
		strncpy(sa.sun_path, name, sizeof(sa.sun_path));
		if (bind(s, (struct sockaddr *) &sa, sizeof(sa)) == -1) {
			perror("ft_start_buffer_server, bind");
			goto cleanup;
		}
		SC->isUnixDomain = 0;
#endif
	} else {
		/* TCP socket */
		struct sockaddr_in sa;

		s = socket(PF_INET, SOCK_STREAM, 0);
		if (s == INVALID_SOCKET) {
			perror("ft_start_buffer_server, socket");
			goto cleanup;
		}
		/* prevend "bind: address already in use" */
		optval = 1;
		if (setsockopt(s, SOL_SOCKET, SO_REUSEADDR, (const char*)&optval, sizeof(optval)) < 0) {
			perror("ft_start_buffer_server, setsockopt");
			goto cleanup;
		}

		bzero(&sa, sizeof(sa));
		sa.sin_family = AF_INET;
		sa.sin_port   = htons(port);
		sa.sin_addr.s_addr = htonl(INADDR_ANY);
		if (bind(s, (struct sockaddr *) &sa, sizeof(sa)) < 0) {
			perror("ft_start_buffer_server, bind");
			goto cleanup;
		}
		SC->isUnixDomain = 0;
	}
	
	/* place the socket in non-blocking mode, required to do thread cancelation */
#ifdef WIN32
	{
		unsigned long enable = 0;
		ioctlsocket(s, FIONBIO, &enable);
	}
#else
	optval = fcntl(s, F_GETFL, NULL);
	optval = optval | O_NONBLOCK;
	if (fcntl(s, F_SETFL, optval)<0) {
		perror("ft_start_buffer_server, fcntl");
		goto cleanup;
	}
#endif

	if (listen(s, BACKLOG)<0) {
		perror("ft_start_buffer_server, listen");
		goto cleanup;
	}
		
	/* set some control variables */
	SC->numClients = 0;
	SC->keepRunning = 1;
	SC->serverSocket = s;
	SC->verbosity = 10; /* TODO: specify proper values */
	
	/* create the mutex */
	if (pthread_mutex_init(&SC->lock, NULL) != 0) {
		fprintf(stderr,"start_tcp_server: mutex could not be initialised\n");
		goto cleanup;
	}
	
	/* create thread with default attributes */
	if (pthread_create(&SC->threadID, NULL, _buffer_server_func, SC) == 0) {
		/* everything went fine - thread should be running now */
		return SC;
	}	
	
	fprintf(stderr,"start_tcp_server: could not spawn thread\n");
	pthread_mutex_destroy(&SC->lock);

cleanup:
	if (SC != NULL) free(SC);
	if (s != INVALID_SOCKET) {
		#ifdef PLATFORM_WIN32
		shutdown(s, SD_BOTH);
		#else
		shutdown(s, SHUT_RDWR);
		#endif
		closesocket(s);
	}
	return NULL;
}		

void ft_stop_buffer_server(ft_buffer_server_t *S) {
	if (S==NULL) return;
	
	S->keepRunning = 0;
	pthread_join(S->threadID, NULL);
	pthread_detach(S->threadID);
	free(S);
}
