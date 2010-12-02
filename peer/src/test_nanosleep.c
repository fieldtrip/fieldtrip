#include <time.h>

void main(void) {
		int retval = 0;
		struct timespec req, rem;
 float t = 1.0001;
		/* split in seconds and nanoseconds */
		req.tv_sec  = (int)t;
		req.tv_nsec = (int)(1000000000.0 * (t - (int)t));
		retval = nanosleep(&req, &rem);
		return retval;
}
