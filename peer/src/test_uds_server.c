#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>

#define SOCK_PATH "echo_socket"

int main(void)
{
		int s, s2, t, len;
		struct sockaddr_un local, remote;
		char str[100];

		if ((s = socket(AF_UNIX, SOCK_STREAM, 0)) == -1) {
				perror("socket");
				exit(1);
		}

		local.sun_family = AF_UNIX;
		strcpy(local.sun_path, SOCK_PATH);
		unlink(local.sun_path);
		len = strlen(local.sun_path) + sizeof(local.sun_family);
		if (bind(s, (struct sockaddr *)&local, len) == -1) {
				perror("bind");
				exit(1);
		}

		if (listen(s, 5) == -1) {
				perror("listen");
				exit(1);
		}

		for(;;) {
				int done, n;
				printf("Waiting for a connection...\n");
				t = sizeof(remote);
				if ((s2 = accept(s, (struct sockaddr *)&remote, &t)) == -1) {
						perror("accept");
						exit(1);
				}

				printf("Connected.\n");

				done = 0;
				do {
						n = recv(s2, str, 100, 0);
						if (n <= 0) {
								if (n < 0) perror("recv");
								done = 1;
						}

						if (!done) 
								if (send(s2, str, n, 0) < 0) {
										perror("send");
										done = 1;
								}
				} while (!done);

				close(s2);
		}

		return 0;
}
