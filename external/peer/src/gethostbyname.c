#include <stdio.h>
#include <stdlib.h>
#include <netdb.h>


int main (void) {
		struct hostent* host;

		while(1) {
				host = gethostbyname("MacBook.local");
		}

		return 0;
}

