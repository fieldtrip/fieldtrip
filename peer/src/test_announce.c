#include <ifaddrs.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <time.h>

#include "peer.h"
#include "extern.h"

int main(int argc, char *argv[]) {

		peerinit(NULL);
		announce(host);
		return 0;
};

