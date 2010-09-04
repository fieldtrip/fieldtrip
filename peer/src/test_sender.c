/*
 * sender.c -- multicasts "hello, world!" to a multicast group once a second
 *
 * Antony Courtney,	25/11/94
 *
 */

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <time.h>
#include <string.h>
#include <stdio.h>

#define HELLO_PORT 12345
#define HELLO_GROUP "225.0.0.37"

main(int argc, char *argv[])
{
		struct sockaddr_in addr;
		int fd, cnt;
		struct ip_mreq mreq;
		char *msgbuf="hello";

		/* create what looks like an ordinary UDP socket */
		if ((fd=socket(AF_INET,SOCK_DGRAM,0)) < 0) {
				perror("socket");
				exit(1);
		}

		/* set up destination address */
		memset(&addr,0,sizeof(addr));
		addr.sin_family=AF_INET;
		addr.sin_addr.s_addr=inet_addr(HELLO_GROUP);
		addr.sin_port=htons(HELLO_PORT);

		/* now just sendto() our destination! */
		while (1) {
				if (sendto(fd,msgbuf,sizeof(msgbuf),0,(struct sockaddr *) &addr,
										sizeof(addr)) < 0) {
						perror("sendto");
						exit(1);
				}
				fprintf(stdout, "%d, %s\n", sizeof(msgbuf), msgbuf);
				sleep(1);
		}
}
