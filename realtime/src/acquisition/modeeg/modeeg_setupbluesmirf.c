/** Set up Bluesmirf (GOLD) baudrate to 57600 temporarily.
	(C) 2010 S. Klanke
*/
#include <stdio.h>
#include <serial.h>

int main(int argc, char *argv[]) {
	SerialPort SP;
	char buf[1025];
	int nr;

	if (argc<1) {
		printf("Usage: modeeg_setupbluesmirf <device>\n");
		printf("  <device> must be your serial port, e.g, COM3: on Windows, or /dev/ttyS0 on Linux\n");
		return 0;
	}
		
	if (!serialOpenByName(&SP, argv[1])) {
		fprintf(stderr, "Could not open serial port %s\n", argv[1]);
		return 1;
	}
	fprintf(stderr, "Opened serial port %s\n", argv[1]);
	
	/* last parameter is timeout in 1/10 of a second */
	if (!serialSetParameters(&SP, 57600, 8, 0, 1, 1)) {
		fprintf(stderr, "Could not modify serial port parameters\n");
		return 1;
	}
	fprintf(stderr, "Modified serial port parameters\n");
	
	nr = serialRead(&SP, 1024, buf);
	
	serialWrite(&SP, 3, "$$$");
	nr = serialRead(&SP, 1024, buf);
	if (nr > 0) {
		buf[nr]=0;
		printf("Received after $$$: %s", buf);
	} else {
		goto cleanup;
	}
	
	serialWrite(&SP, 2, "D\r");
	nr = serialRead(&SP, 1024, buf);
	if (nr > 0) {
		buf[nr]=0;
		printf("Received after D<cr>: %s", buf);
	} else {
		goto cleanup;
	}
	
	serialWrite(&SP, 9, "U,576K,N\r");
	nr = serialRead(&SP, 1024, buf);
	if (nr > 0) {
		buf[nr]=0;
		printf("Received after U,576K,N<cr>: %s", buf);
	} else {
		goto cleanup;
	}	
	
cleanup:
	serialClose(&SP);
	
	return 0;
}
