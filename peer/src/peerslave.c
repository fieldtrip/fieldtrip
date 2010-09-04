/*
 * Copyright (C) 2010, Robert Oostenveld
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/
 *
 */ 

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>

#include "engine.h"

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

mxArray *mxSerialize(const mxArray*);
mxArray *mxDeserialize(const void*, size_t);

#define ENGINETIMEOUT     30     /* in seconds */
#define ZOMBIETIMEOUT     300    /* in seconds */
#define SLEEPTIME         10000  /* in microseconds */

#define STARTCMD "matlab -nosplash"

void print_help(char *argv[]) {
		printf("\n");
		printf("This starts a FieldTrip peer-to-peer distributed computing peer, which\n");
		printf("will wait for an incoming job and subsequently start the MATLAB engine and\n");
		printf("evaluate the job. Use as\n");
		printf("  %s [options]\n", argv[0]);
		printf("where the options can include\n");
		printf("  --memavail    = number, amount of memory available       (default = inf)\n");
		printf("  --cpuavail    = number, speed of the CPU                 (default = inf)\n");
		printf("  --timavail    = number, maximum duration of a single job (default = inf)\n");
		printf("  --allowhost   = {...}\n");
		printf("  --allowuser   = {...}\n");
		printf("  --allowgroup  = {...}\n");
		printf("  --group       = string\n");
		printf("  --hostname    = string\n");
		printf("  --matlab      = string\n");
		printf("  --timeout     = number, time to keep the engine running after the job finished\n");
		printf("  --smartshare  = 0|1\n");
		printf("  --smartmem    = 0|1\n");
		printf("  --smartcpu    = 0|1\n");
		printf("  --daemon\n");
		printf("  --udsserver\n");
		printf("  --verbose\n");
		printf("  --help\n");
		printf("\n");
}

int help_flag;
int daemon_flag;
int verbose_flag;
int tcpserver_flag = 1;
int udpserver_flag = 0;
int udsserver_flag = 0;

int main(int argc, char *argv[]) {
		Engine *en;
		mxArray *argin, *argout, *options, *arg, *opt;
		joblist_t  *job  = NULL;
		peerlist_t *peer = NULL;
		jobdef_t   *def  = NULL;
		pid_t childpid;

		int matlabRunning = 0, matlabStart, matlabFinished, engineFailed = 0;
		int c, n, rc, found, handshake, success, server, jobnum = 0;
		unsigned int enginetimeout = ENGINETIMEOUT;
		unsigned int zombietimeout = ZOMBIETIMEOUT;
		unsigned int peerid, jobid;
		char *str = NULL, *startcmd = NULL;

		userlist_t  *allowuser;
		grouplist_t *allowgroup;
		hostlist_t  *allowhost;

		/* the thread IDs are needed for cancelation at cleanup */
		pthread_t udsserverThread;
		pthread_t tcpserverThread;
		pthread_t announceThread;
		pthread_t discoverThread;
		pthread_t expireThread;

		openlog("peerslave", 0, LOG_USER);
		peerinit(NULL);

		/* use GNU getopt_long for the command-line options */
		while (1)
		{
				static struct option long_options[] =
				{
						{"help",       no_argument, &help_flag, 1},
						{"daemon",     no_argument, &daemon_flag, 1},
						{"verbose",    no_argument, &verbose_flag, 1},
						{"udsserver",  no_argument, &udsserver_flag, 1},
						{"memavail",   required_argument, 0, 'a'}, /* numeric argument */
						{"cpuavail",   required_argument, 0, 'b'}, /* numeric argument */
						{"timavail",   required_argument, 0, 'c'}, /* numeric argument */
						{"hostname",   required_argument, 0, 'd'}, /* single string argument */
						{"group",      required_argument, 0, 'e'}, /* single string argument */
						{"allowuser",  required_argument, 0, 'f'}, /* single or multiple string argument */
						{"allowhost",  required_argument, 0, 'g'}, /* single or multiple string argument */
						{"allowgroup", required_argument, 0, 'h'}, /* single or multiple string argument */
						{"matlab",     required_argument, 0, 'i'}, /* single string argument */
						{"smartmem",   required_argument, 0, 'j'}, /* numeric, 0 or 1 */
						{"smartcpu",   required_argument, 0, 'k'}, /* numeric, 0 or 1 */
						{"smartshare",  required_argument, 0, 'l'}, /* numeric, 0 or 1 */
						{"timeout",    required_argument, 0, 'm'}, /* numeric argument */
						{0, 0, 0, 0}
				};

				/* getopt_long stores the option index here. */
				int option_index = 0;

				c = getopt_long (argc, argv, "", long_options, &option_index); 

				/* detect the end of the options. */
				if (c == -1)
						break;

				switch (c)
				{
						case 0:
								/* if this option set a flag, do nothing else now. */
								if (long_options[option_index].flag != 0)
										break;
								break;

						case 'a':
								if (verbose_flag)
										printf ("option --memavail with value `%s'\n", optarg);
								pthread_mutex_lock(&mutexhost);
								host->memavail = atol(optarg);
								pthread_mutex_unlock(&mutexhost);
								pthread_mutex_lock(&mutexsmartmem);
								smartmem.enabled = 0;
								pthread_mutex_unlock(&mutexsmartmem);
								break;

						case 'b':
								if (verbose_flag)
										printf ("option --cpuavail with value `%s'\n", optarg);
								pthread_mutex_lock(&mutexhost);
								host->cpuavail = atol(optarg);
								pthread_mutex_unlock(&mutexhost);
								break;

						case 'c':
								if (verbose_flag)
										printf ("option --timavail with value `%s'\n", optarg);
								pthread_mutex_lock(&mutexhost);
								host->timavail = atol(optarg);
								pthread_mutex_unlock(&mutexhost);
								break;

						case 'd':
								if (verbose_flag)
										printf ("option --hostname with value `%s'\n", optarg);
								pthread_mutex_lock(&mutexhost);
								strncpy(host->name, optarg, STRLEN);
								pthread_mutex_unlock(&mutexhost);
								break;

						case 'e':
								if (verbose_flag)
										printf ("option --group with value `%s'\n", optarg);
								pthread_mutex_lock(&mutexhost);
								strncpy(host->group, optarg, STRLEN);
								pthread_mutex_unlock(&mutexhost);
								break;

						case 'f':
								if (verbose_flag)
										printf ("option --allowuser with value `%s'\n", optarg);
								str = strtok(optarg, ",");
								while (str) {
										allowuser = (userlist_t *)malloc(sizeof(userlist_t));
										allowuser->name = malloc(STRLEN);
										strncpy(allowuser->name, str, STRLEN);
										allowuser->next = userlist;
										userlist = allowuser;
										str = strtok(NULL, ",");
								}
								break;

						case 'g':
								if (verbose_flag)
										printf ("option --allowhost with value `%s'\n", optarg);
								str = strtok(optarg, ",");
								while (str) {
										allowhost = (hostlist_t *)malloc(sizeof(hostlist_t));
										allowhost->name = malloc(STRLEN);
										strncpy(allowhost->name, str, STRLEN);
										allowhost->next = hostlist;
										hostlist = allowhost;
										str = strtok(NULL, ",");
								}
								break;

						case 'h':
								if (verbose_flag)
										printf ("option --allowgroup with value `%s'\n", optarg);
								str = strtok(optarg, ",");
								while (str) {
										allowgroup = (grouplist_t *)malloc(sizeof(grouplist_t));
										allowgroup->name = malloc(STRLEN);
										strncpy(allowgroup->name, str, STRLEN);
										allowgroup->next = grouplist;
										grouplist = allowgroup;
										str = strtok(NULL, ",");
								}
								break;

						case 'i':
								if (verbose_flag)
										printf ("option --matlab with value `%s'\n", optarg);
								startcmd = malloc(STRLEN);
								strncpy(startcmd, optarg, STRLEN);
								break;

						case 'j':
								if (verbose_flag)
										printf ("option --smartmem with value `%s'\n", optarg);
								pthread_mutex_lock(&mutexsmartmem);
								smartmem.enabled = atol(optarg);
								pthread_mutex_unlock(&mutexsmartmem);
								break;

						case 'k':
								if (verbose_flag)
										printf ("option --smartcpu with value `%s'\n", optarg);
								pthread_mutex_lock(&mutexsmartcpu);
								smartcpu.enabled = atol(optarg);
								pthread_mutex_unlock(&mutexsmartcpu);
								break;

						case 'l':
								if (verbose_flag)
										printf ("option --smartshare with value `%s'\n", optarg);
								pthread_mutex_lock(&mutexsmartshare);
								smartshare.enabled = atol(optarg);
								pthread_mutex_unlock(&mutexsmartshare);

						case 'm':
								if (verbose_flag)
										printf ("option --timeout with value `%s'\n", optarg);
								enginetimeout = atol(optarg);
								break;

								break;

						case '?':
								/* getopt_long already printed an error message. */
								break;

						default:
								fprintf(stderr, "invalid command line options\n");
								exit(1);
				}
		}

		if (help_flag) {
				/* display the help message and return to the command line */
				print_help(argv);
				exit(0);
		}

		if (daemon_flag) {
				/* now create new process */
				childpid = fork();

				if (childpid >= 0) /* fork succeeded */
				{
						if (childpid == 0) /* fork() returns 0 to the child process */
						{
								printf("CHILD: I am the child process!\n");
								printf("CHILD: Here's my PID: %d\n", getpid());
								printf("CHILD: My parent's PID is: %d\n", getppid());
								printf("CHILD: The value of my copy of childpid is: %d\n", childpid);
								/* the child continues as the actual executable */
						}
						else /* fork() returns new pid to the parent process */
						{
								printf("PARENT: I am the parent process!\n");
								printf("PARENT: Here's my PID: %d\n", getpid());
								printf("PARENT: The value of my copy of childpid is %d\n", childpid);
								/* parent exits */ 
								return 0;
						}
				}
				else /* fork returns -1 on failure */
				{
						perror("fork"); /* display error message */
						syslog(LOG_ERR, "error: fork"); /* display error message */
						exit(0); 
				}
		}

		/* set the default command to start matlab */
		if (!startcmd) {
				startcmd = malloc(STRLEN);
				strncpy(startcmd, STARTCMD, STRLEN);
		}

		if (udsserver_flag) {
				if ((rc = pthread_create(&udsserverThread, NULL, udsserver, (void *)NULL))>0) {
						fprintf(stderr, "failed to start udsserver thread\n");
						exit(1);
				}
				else {
						if (verbose_flag)
								fprintf(stderr, "started udsserver thread\n");
				}
		}

		if (tcpserver_flag) {
				if ((rc = pthread_create(&tcpserverThread, NULL, tcpserver, (void *)NULL))>0) {
						fprintf(stderr, "failed to start tcpserver thread\n");
						exit(1);
				}
				else {
						if (verbose_flag)
								fprintf(stderr, "started tcpserver thread\n");
				}
		}

		if ((rc = pthread_create(&announceThread, NULL, announce, (void *)NULL))>0) {
				fprintf(stderr, "failed to start announce thread\n");
				exit(1);
		}
		else {
				if (verbose_flag)
						fprintf(stderr, "started announce thread\n");
		}

		if ((rc = pthread_create(&discoverThread, NULL, discover, (void *)NULL))>0) {
				fprintf(stderr, "failed to start discover thread\n");
				exit(1);
		}
		else {
				if (verbose_flag)
						fprintf(stderr, "started discover thread\n");
		}

		if ((rc = pthread_create(&expireThread, NULL, expire, (void *)NULL))>0) {
				fprintf(stderr, "failed to start expire thread\n");
				exit(1);
		}
		else {
				if (verbose_flag)
						fprintf(stderr, "started expire thread\n");
		}

		/* status = 0 means zombie mode, don't accept anything   */
		/* status = 1 means master mode, accept everything       */
		/* status = 2 means idle slave, accept only a single job */
		/* status = 3 means busy slave, don't accept a new job   */
		/* any other status is interpreted as zombie mode        */
		host->status = 2;

		while (1) {

				if (jobcount()>0 && !engineFailed) {

						/* there is a job to be executed */
						if (matlabRunning==0) {
								fprintf(stderr, "starting MATLAB engine\n");
								if ((en = engOpen(startcmd)) == NULL) {
										/* this may be due to a licensing problem */
										/* do not attempt to start again during the timeout period */
										pthread_mutex_lock(&mutexhost);
										fprintf(stderr, "failed to start MATLAB engine\n");
										fprintf(stderr, "switching to zombie mode\n");
										engineFailed = time(NULL);
										host->status = 0; /* zombie */
										pthread_mutex_unlock(&mutexhost);

										/* remove the first job from the joblist */
										fprintf(stderr, "deleting job\n");
										pthread_mutex_lock(&mutexjoblist);
										job = joblist;
										joblist = job->next; 
										FREE(job->job);
										FREE(job->host);
										FREE(job->arg);
										FREE(job->opt);
										FREE(job);
										pthread_mutex_unlock(&mutexjoblist);
										continue;
								}
								else {
										matlabRunning = 1;
								}
						}

						/* switch the mode to busy slave */
						host->status = 3; /* busy */
						matlabStart = time(NULL);

						/* get the first job input arguments and options */
						pthread_mutex_lock(&mutexjoblist);
						job = joblist;
						argin   = (mxArray *)mxDeserialize(job->arg, job->job->argsize);
						options = (mxArray *)mxDeserialize(job->opt, job->job->optsize);
						jobid   = job->job->id;
						peerid  = job->host->id;
						fprintf(stderr, "executing job %d from %s@%s (jobid=%d)\n", ++jobnum, job->host->user, job->host->name, jobid);
						pthread_mutex_unlock(&mutexjoblist);

						/* copy them over to the engine */
						engPutVariable(en, "argin", argin);
						engPutVariable(en, "options", options);
						mxDestroyArray(argin);
						mxDestroyArray(options);

						/* execute the job */
						engEvalString(en, "[argout, options] = peerexec(argin, options);");

						/* get the job output arguments and options */
						if ((argout = engGetVariable(en, "argout")) == NULL) {
								fprintf(stderr, "error getting argout");
								exit(1);
						}

						if ((options = engGetVariable(en, "options")) == NULL) {
								fprintf(stderr, "error getting options");
								exit(1);
						}

						/* send them back to the master */

						/*****************************************************************************
						 * the following code is largely shared with the get-option in the peer mex file
						 *****************************************************************************/
						pthread_mutex_lock(&mutexpeerlist);
						found = 0;

						peer = peerlist;
						while(peer) {
								if (peer->host->id==peerid) {
										found = 1;
										break;
								}
								peer = peer->next ;
						}

						if (!found) {
								pthread_mutex_unlock(&mutexpeerlist);
								PANIC("failed to locate specified peer\n");
						}

						pthread_mutex_lock(&mutexhost);
						int hasuds = (strlen(peer->host->socket)>0 && strcmp(peer->host->name, host->name)==0);
						int hastcp = (peer->host->port>0);
						pthread_mutex_unlock(&mutexhost);

						if (hasuds) {
								/* open the UDS socket */
								if ((server = open_uds_connection(peer->host->socket)) < 0) {
										pthread_mutex_unlock(&mutexpeerlist);
										PANIC("failed to create socket\n");
								}
						}
						else if (hastcp) {
								/* open the TCP socket */
								if ((server = open_tcp_connection(peer->ipaddr, peer->host->port)) < 0) {
										pthread_mutex_unlock(&mutexpeerlist);
										PANIC("failed to create socket\n");
								}
						}
						else {
								pthread_mutex_unlock(&mutexpeerlist);
								PANIC("failed to create socket\n");
						}

						/* the connection was opened without error */
						pthread_mutex_unlock(&mutexpeerlist);

						if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
								PANIC("could not write handshake");
						}
						else if (!handshake) {
								close_connection(server);
								PANIC("failed to negociate connection");
						}

						/* the message that will be written consists of
						   message->host
						   message->job
						   message->arg
						   message->opt
						 */

						arg = (mxArray *) mxSerialize(argout);
						opt = (mxArray *) mxSerialize(options);

						if (!arg)
								PANIC("could not serialize job arguments");

						if (!opt)
								PANIC("could not serialize job options");

						def = (jobdef_t *)malloc(sizeof(jobdef_t));

						if (!def)
								PANIC("could not allocate memory");

						def->version  = VERSION;
						def->id       = jobid;
						def->memreq   = 0;
						def->cpureq   = 0;
						def->timreq   = 0;
						def->argsize  = mxGetNumberOfElements(arg);
						def->optsize  = mxGetNumberOfElements(opt);

						/* write the message  (hostdef, jobdef, arg, opt) with handshakes in between */
						/* the slave may close the connection between the message segments in case the job is refused */
						success = 1;

						pthread_mutex_lock(&mutexhost);
						if (success) 
								success = (bufwrite(server, host, sizeof(hostdef_t)) == sizeof(hostdef_t));
						pthread_mutex_unlock(&mutexhost);

						if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
								PANIC("could not write handshake");
						}
						else if (!handshake) {
								close_connection(server);
								PANIC("failed to write hostdef");
						}

						if (success)
								success = (bufwrite(server, def, sizeof(jobdef_t)) == sizeof(jobdef_t));

						if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
								PANIC("could not write handshake");
						}
						else if (!handshake) {
								close_connection(server);
								PANIC("failed to write jobdef");
						}

						if (success) 
								success = (bufwrite(server, (void *)mxGetData(arg), def->argsize) == def->argsize);

						if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
								PANIC("could not write handshake");
						}
						else if (!handshake) {
								close_connection(server);
								PANIC("failed to write arg");
						}

						if (success) 
								success = (bufwrite(server, (void *)mxGetData(opt), def->optsize) == def->optsize);

						if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
								PANIC("could not write handshake");
						}
						else if (!handshake) {
								close_connection(server);
								PANIC("failed to write opt");
						}

						close_connection(server);

						mxDestroyArray(arg);
						mxDestroyArray(opt);

						/*****************************************************************************/

						/* remove the first job from the joblist */
						pthread_mutex_lock(&mutexjoblist);
						joblist = job->next;
						FREE(job->job);
						FREE(job->host);
						FREE(job->arg);
						FREE(job->opt);
						FREE(job);
						job = joblist;
						pthread_mutex_unlock(&mutexjoblist);

						host->status = 2; /* idle slave */
						matlabFinished = time(NULL);
						fprintf(stderr, "executing job %d took %d seconds\n", jobnum, matlabFinished - matlabStart);
				}
				else {
						usleep(SLEEPTIME);
				} /* if jobcount */

				/* test that the matlab engine is not idle for too long */
				if ((matlabRunning==1) && (time(NULL)-matlabFinished)>enginetimeout) {
						if (engClose(en)!=0) {
								PANIC("could not stop the MATLAB engine\n");
						}
						else {
								fprintf(stderr, "stopped idle MATLAB engine\n");
								matlabRunning = 0;
						}
				}

				/* don't try to restart the engine immediately after a failure */
				if (engineFailed && (time(NULL)-engineFailed)>zombietimeout) {
						pthread_mutex_lock(&mutexhost);
						fprintf(stderr, "switching to idle mode\n");
						engineFailed = 0;
						host->status = 2;  /* idle slave */
						pthread_mutex_unlock(&mutexhost);
						sleep(1);
						continue;
				}

		} /* while */

		return 0;
} /* main */

