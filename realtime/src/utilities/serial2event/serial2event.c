/* Convert incoming characters on a serial port to an event in a remote FieldTrip buffer.
 * Errors when writing to the buffer will be printed, but otherwise ignored.
 * Please look at serial2event.conf for an example of how to set up the tool,
 * e.g., which port to listen on, and how to write events.
 * (C) 2010 Stefan Klanke
 */
#include <serial.h>
#include <buffer.h>
#include <signal.h>
#include <pthread.h>

#define MAXLINE 256

typedef struct {
	char hostname[256];
	int port;
	
	int character; /* 0..127 for a certain character, or -1 for don't care */
	int set_value; /* non-zero: transmit incoming character as "value" field */
	
	char comport[256];
	int baudrate, databits, stopbits, parity;
	
	int udp_port;
	
	INT32_T sample_start;
	INT32_T sample_increase;
	INT32_T offset, duration;
	UINT32_T type_type;
	UINT32_T type_numel;
	UINT32_T value_type;
	UINT32_T value_numel;
	
	char type_buf[MAXLINE];
	char value_buf[MAXLINE];
} SerialEventConfig;



int keepRunning = 1;
int udp_socket;
SerialEventConfig conf;
INT32_T sample;
pthread_mutex_t sampleMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_t udpThread;


/* The following function is used for turning a configuration file (see serial2event.conf)
   into a SerialEventConfig struct as defined above.
   Returns 
	 0 if ok, 
	-1 if C==NULL, 
	-2 if file can't be opened, 
	positive number = #errors 
*/
int parseConfig(SerialEventConfig *C, const char *filename) {
	FILE *f;
	char line[MAXLINE];
	int numErrs = 0;
	int lineNr = 0;
	
	if (C==NULL) return -1;
	
	/* set defaults */
	strcpy(C->hostname, "localhost");
	C->port = 1972;
	C->character = -1;
	C->set_value = 1;
	C->sample_start = 0;
	C->sample_increase = 1;
	C->offset = C->duration = 0;
	C->type_type = DATATYPE_CHAR;
	C->type_numel = 6;
	strcpy(C->type_buf, "serial");
	C->value_type = DATATYPE_CHAR;
	C->type_numel = 1;
	C->udp_port = 1990;
	
	f = fopen(filename, "r");
	if (f==NULL) {
		printf("Configuration file %s could not be opened\n", filename);
		return -2;
	}
	
	while (!feof(f)) {
		int len;
		char *value;
		
		if (fgets(line, MAXLINE, f) == NULL) break;
		
		++lineNr;
		
		/* silently ignore comments */
		if (line[0]=='#') continue;
		
		len = strlen(line);
		/* strip trailing newline, carriage return */
		while (len>0) {
			int c = line[len-1];
			if (c=='\r' || c=='\n') {
				--len;
			} else {
				break;
			}
		}
		line[len]=0;
		
		/* silently ignore empty lines */
		if (len==0) continue;
		
		value = strchr(line,'=');
		if (value == NULL) {
			++numErrs;
			printf("Ignoring faulty line %i\n", lineNr);
			continue;
		}
		*value++ = 0;
		if (!strcmp(line, "buffer")) {
			char *pc;
			long lv = 0;
			
			pc = strchr(value, ':');
			if (pc != NULL) {
				lv = strtol(pc+1, NULL, 10);
			}
			if (lv<=0 || pc == value) {
				printf("Ignoring faulty buffer target defintion at line %i\n", lineNr);
				++numErrs;
			} else {
				*pc = '\0'; /* replace ':' by terminator */
				strcpy(C->hostname, value);
				C->port = lv;
			}
		} else if (!strcmp(line, "character")) {
			C->character = *value;
			printf("Filtering for serial port character [%c]\n", C->character);
		} else if (!strcmp(line, "sample")) {
			char *p1, *p2;
			long lva,lvb;
	
			p1 = strchr(value, '+');
			if (p1!=NULL) {
				*p1++ = 0; /* replace '+' by '\0' */
			}

			lva = strtol(value, &p2, 10);
			if (p2==value) {
				++numErrs;
				printf("Ignoring faulty sample index in line %i\n", lineNr);
				continue;
			}
			
			if (p1==NULL) {
				C->sample_increase = 0;
				C->sample_start = (INT32_T) lva;
				continue;
			}
			
			lvb = strtol(p1, &p2, 10);
			if (p2==p1) {
				++numErrs;
				printf("Ignoring faulty sample increase in line %i\n", lineNr);
				continue;
			}
			C->sample_start = lva;
			C->sample_increase = lvb;
		} else if (!strcmp(line, "port")) {
			char *p;
			long lv;
		
			lv = strtol(value, &p, 10);
			if (value==p) {
				++numErrs;
				printf("Ignoring faulty 'port' in line %i\n", lineNr);
				continue;
			}
			C->udp_port = lv;
		} else if (!strcmp(line, "offset")) {
			char *p;
			long lv;
		
			lv = strtol(value, &p, 10);
			if (value==p) {
				++numErrs;
				printf("Ignoring faulty offset in line %i\n", lineNr);
				continue;
			}
			C->offset = lv;
		} else if (!strcmp(line, "duration")) {
			char *p;
			long lv;
		
			lv = strtol(value, &p, 10);
			if (value==p) {
				++numErrs;
				printf("Ignoring faulty duration in line %i\n", lineNr);
				continue;
			}
			C->duration = lv;
		} else if (!strcmp(line, "serial")) {
			int baud, data, stop, par, num;
			char *p;
			
			p = strchr(value, ':');
			if (p==NULL) {
				num = 0;
			} else {
				*p++ = 0; /* replace first ':' by terminator */
				num = sscanf(p, "%i:%i:%i:%i", &baud, &data, &stop, &par);
			}
			if (num<4 || par<0 || par>1 || stop<0 || data<0 || data>9 || stop>2 || baud<0) {
				printf("Ignoring faulty serial parameter defintion in line %i\n", lineNr);
				++numErrs;
			} else {
				strcpy(C->comport, value);
				C->baudrate = baud;
				C->databits = data;
				C->stopbits = stop;
				C->parity = par;
			}
		} else if (!strcmp(line, "type")) {
			char *p1,*p2;
			p1 = strchr(value, '"');
			if (p1!=NULL) {
				p1++;
				p2 = strchr(p1, '"');
				if (p2==NULL) {
					++numErrs;
					printf("Unterminated string in type field, line %i\n", lineNr);
				} else {
					C->type_numel = p2-p1;
					C->type_type = DATATYPE_CHAR;
					memcpy(C->type_buf, p1, C->type_numel);
				}
			} else {
				p1 = strchr(value, '.');
				if (p1!=NULL) {
					/* double value */
					double dv;
					dv = strtod(value, &p2);
					if (p2==value) {
						printf("Invalid double precision value for type field, line %i\n", lineNr);
						++numErrs;
					} else {
						C->type_numel = 1;
						C->type_type  = DATATYPE_FLOAT64;
						memcpy(C->type_buf, &dv, sizeof(double));
					}
				} else {
					/* int value */
					INT32_T iv;
					iv = (INT32_T) strtol(value, &p2, 10);
					if (p2==value) {
						printf("Invalid integer value for type field, line %i\n", lineNr);
						++numErrs;
					} else {
						C->type_numel = 1;
						C->type_type  = DATATYPE_INT32;
						memcpy(C->type_buf, &iv, sizeof(iv));
					}
				}
			}
		} else if (!strcmp(line, "value")) {
			if (*value == '@') {
				C->value_type = DATATYPE_CHAR;
				C->value_numel = 1;
				C->value_buf[0] = '@';
				C->set_value = 1;
			} else {
				char *p1,*p2;
				p1 = strchr(value, '"');
				if (p1!=NULL) {
					p1++;
					p2 = strchr(p1, '"');
					if (p2==NULL) {
						++numErrs;
						printf("Unterminated string in value field, line %i\n", lineNr);
					} else {
						C->value_type = DATATYPE_CHAR;
						C->value_numel = p2-p1;
						memcpy(C->value_buf, p1, C->value_numel);
						C->set_value = 0;
					}
				} else {
					p1 = strchr(value, '.');
					if (p1!=NULL) {
						/* double value */
						double dv;
						dv = strtod(value, &p2);
						if (p2==value) {
							printf("Invalid double precision value for type field, line %i\n", lineNr);
							++numErrs;
						} else {
							C->value_numel = 1;
							C->value_type  = DATATYPE_FLOAT64;
							memcpy(C->value_buf, &dv, sizeof(double));
							C->set_value = 0;
						}
					} else {
						/* int value */
						INT32_T iv;
						iv = (INT32_T) strtol(value, &p2, 10);
						if (p2==value) {
							printf("Invalid integer value for type field, line %i\n", lineNr);
							++numErrs;
						} else {
							C->value_numel = 1;
							C->value_type  = DATATYPE_INT32;
							memcpy(C->value_buf, &iv, sizeof(iv));
							C->set_value = 0;
						}
					}
				}
			}
		
		} else {
			printf("Ignoring unknown field defintion in line %i\n", lineNr);
		}
	}
	fclose(f);
	
	printf("\nEvent definition\n----------------\n");
	printf("Sample....: %i + %i\n", C->sample_start, C->sample_increase);
	printf("Offset....: %i\n", C->offset);
	printf("Duration..: %i\n", C->duration);
	printf("Type......: ");
	switch(C->type_type) {
		case DATATYPE_CHAR:
			printf("'%.*s' (%i characters)\n", C->type_numel, C->type_buf, C->type_numel);
			break;
		case DATATYPE_INT32:
			printf("%i (int32)\n", *((INT32_T *) C->type_buf));
			break;
		case DATATYPE_FLOAT64:
			printf("%f (float64)\n", *((double *) C->type_buf));
			break;
	}
	printf("Value.....: ");
	if (C->set_value) {
		printf("<from serial port>\n");
	} else switch(C->value_type) {
		case DATATYPE_CHAR:
			printf("'%.*s' (%i characters)\n", C->value_numel, C->value_buf, C->value_numel);
			break;
		case DATATYPE_INT32:
			printf("%i (int32)\n", *((INT32_T *) C->value_buf));
			break;
		case DATATYPE_FLOAT64:
			printf("%f (float64)\n", *((double *) C->value_buf));
			break;
	}
	
	printf("\nSerial port\n-----------\n");
	printf("Device....: %s\n", C->comport);
	printf("Baudrate..: %i\n", C->baudrate);
	printf("Data bits.: %i\n", C->databits);
	printf("Stop bits.: %i\n", C->stopbits);
	printf("Parity....: %i\n", C->parity);
	if (C->character == -1) {
		printf("No filter\n");
	} else {
		printf("React on..: %c\n", C->character);
	}
	printf("\nWill write to buffer at %s:%i\n\n", C->hostname, C->port);
	return numErrs;
}


/* Creates a UDP socket and starts listening on the given port number.
   The purpose of this is to detect "RESET" messages (5 characters)
   that might be send, e.g., from the gui_streamer fMRI export tool,
   such that this tool knows that the internal sample counter should
   be reset.
*/
int create_udp_receiver(int port) {
	struct sockaddr_in addr;
	int sock;
	unsigned long	optval;

	memset(&addr, 0, sizeof(struct sockaddr_in));
	addr.sin_family = AF_INET;
	addr.sin_port = htons(port);
	addr.sin_addr.s_addr = htonl(INADDR_ANY);

	sock = socket(PF_INET, SOCK_DGRAM, IPPROTO_UDP);
	if (sock == -1) {
		perror("socket(...): ");
		return -1;
	}
	if (bind(sock, (struct sockaddr *) &addr, sizeof(struct sockaddr_in))) {
		perror("bind(...): ");
		closesocket(sock);
		return -1;
	}
	
#ifdef PLATFORM_WIN32
	optval = 1;
	if (ioctlsocket(sock, FIONBIO, &optval) != 0) {
		fprintf(stderr,"Could not set non-blocking mode\n");
		closesocket(sock);
		return -1;
	}
#else
	optval = fcntl(sock, F_GETFL, NULL);
	optval = optval | O_NONBLOCK;
	if (fcntl(sock, F_SETFL, optval)<0) {
		perror("fcntl (set non-blocking)");
		closesocket(sock);
		return -1;
	}
#endif	
	return sock;
}

/** Main function of background thread that listens for UDP messages
    on the socket that has been setup by the function above. Resets
	the sample counter to conf.sample_start in case a "RESET" message
	comes in. Note that there is no further synchronisation between
	this thread and the main thread, that is, "sample" is just overwritten
	here. This is ok because the operation is atomic.
*/
void *_udp_thread(void *arg) {
	fd_set readSet;
	struct timeval tv;
	char msg[8192];
	
	tv.tv_sec  = 0; 
	tv.tv_usec = 100000; /* 100 ms  timeout for select */
	
	while (keepRunning) {
		int n;
		
		FD_ZERO(&readSet);
		FD_SET(udp_socket, &readSet);
		
		n = select(udp_socket + 1, &readSet, NULL, NULL, &tv);
		if (n<1) continue;
		
		n = recv(udp_socket, msg, sizeof(msg)-1, 0);
		if (n<=0) continue;
		
		/* force-terminate string */
		msg[n] = 0;
		
		printf("Message received: %s\n", msg);
		if (!strcmp(msg,"RESET")) {
			pthread_mutex_lock(&sampleMutex);
			sample = conf.sample_start;
			pthread_mutex_unlock(&sampleMutex);
		}	
	}
	return NULL;
}


/** Function that is called when the user presses CTRL-C */
void abortHandler(int sig) {
	printf("Ctrl-C pressed -- stopping...\n");
	keepRunning = 0;
}

int main(int argc, char **argv) {
	SerialPort SP;
	int ftBuffer = -1;
	eventdef_t *evdef;
	UINT32_T sizetype, sizevalue, bufsize;
	char *valBuf;
	messagedef_t reqdef;
	message_t request, *response;
	char *confname;
	
	if (argc < 2) {
		confname = "serial2event.conf";
	} else {
		confname = argv[1];
	}
	
	if (parseConfig(&conf, confname) != 0) {
		printf("Errors during parsing the configuration file\n");
		exit(1);
	}
	
	sizetype = wordsize_from_type(conf.type_type) * conf.type_numel;
	sizevalue = wordsize_from_type(conf.value_type) * conf.value_numel;
	
	bufsize = sizeof(eventdef_t) + sizetype + sizevalue;
	evdef = (eventdef_t *) malloc(bufsize);
	
	if (evdef == NULL) {
		printf("Out of memory\n");
		exit(1);
	}
	
	/* prepare fixed fields */
	reqdef.version = VERSION;
	reqdef.bufsize = bufsize;
	reqdef.command = PUT_EVT;
	request.def = &reqdef;
	request.buf = evdef;
	
	evdef->offset = conf.offset;
	evdef->duration = conf.duration;
	evdef->type_type = conf.type_type;
	evdef->type_numel = conf.type_numel;
	evdef->value_type = conf.value_type;
	evdef->value_numel = conf.value_numel;
	evdef->bufsize = sizetype + sizevalue;
	
	valBuf = (char *) evdef + sizeof(eventdef_t);
	memcpy(valBuf, conf.type_buf, sizetype);
	valBuf += sizetype;
	memcpy(valBuf, conf.value_buf, sizevalue);
	
	if (!serialOpenByName(&SP, conf.comport)) {
		printf("Could not open serial port.\n");
		exit(1);
	}
	
	/* timeout = 1 decisecond (least common denominator) = 100 ms */
	if (!serialSetParameters(&SP, conf.baudrate, conf.databits, conf.parity, conf.stopbits, 1)) {
		printf("Could not modify serial port parameters\n");
		exit(1);
	}
	
	ftBuffer = open_connection(conf.hostname, conf.port);
	if (ftBuffer < 0) {
		printf("Connection to FieldTrip buffer failed.\n");
		exit(1);
	}
	
	sample = conf.sample_start;
	
	udp_socket = create_udp_receiver(conf.udp_port);
	
	if (udp_socket != -1) {
		if (pthread_create(&udpThread, NULL, _udp_thread, NULL)) {
			printf("Warning: UDP socket thread could not be spawned.\n");
			closesocket(udp_socket);
			udp_socket = -1;
		}
	}
	
	/* register CTRL-C handler */
	signal(SIGINT, abortHandler);
	
	printf("Starting to listen - press CTRL-C to quit\n");
	while (keepRunning) {
		char input;
		int n;
				
		n = serialRead(&SP, 1, &input);
		if (n<0) {
			printf("Error while reading from serial port - exiting\n");
			break;
		}
		if (n==0) continue; /* timeout - just wait longer */
		
		/* we got one */
		if (conf.character != -1 && conf.character != input) {
			printf("Ignoring input %c\n", input);
			continue;
		}
		
		if (conf.set_value) *valBuf = input;
		pthread_mutex_lock(&sampleMutex);
		evdef->sample = sample;
		pthread_mutex_unlock(&sampleMutex);
		
		if (evdef->sample < 0) {
			printf("Ignoring negative sample (%i) event...\n", evdef->sample);
		} else {
			n = tcprequest(ftBuffer, &request, &response);
		
			if (n<0 || response == NULL) {
				printf("Error in FieldTrip connection\n");
			} else {
				if (response->def == NULL || response->def->command != PUT_OK) {
					printf("FieldTrip server returned an error\n");
				} else {
					printf("Sent off event (sample = %i, input = %c)\n", evdef->sample, input);
				}
				FREE(response->def);
				FREE(response->buf);
				free(response);
			}
		}
		
		pthread_mutex_lock(&sampleMutex);
		sample += conf.sample_increase;
		pthread_mutex_unlock(&sampleMutex);
	}
	
	if (udp_socket!=-1) {
		pthread_join(udpThread, NULL);
		closesocket(udp_socket);
	}

	close_connection(ftBuffer);
	serialClose(&SP);
	free(evdef);
	
	return 0;
}
