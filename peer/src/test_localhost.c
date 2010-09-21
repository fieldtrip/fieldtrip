#include <stdio.h>
#include "peer.h"

int main(int argc, char *argv[]) {
		int val;

		while (argc>1) {
				val = check_localhost(argv[1]);
				printf("val = %d\n", val);
		}

		return 0;

}
