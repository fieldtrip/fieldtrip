/** Small class for getting a measurement of the system clock in seconds
    and fractions thereof. On Windows, Multimedia timers are used to get
    millisecond precision, thus you need to link against winmm.lib (or
	libwinmm.a for MinGW). Other platforms are handled with the
	standard 'gettimeofday'.
*/

#ifndef __Clock_h
#define __Clock_h

#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

class Clock {
	public:

	/** Constructor. Will set internal time base to current time, so
		that succesive getRel() calls will return the time lapsed since
		construction of this object.
		On Windows, this initialises the multimedia timer to millisecond
		precision.
	*/
	Clock() {
		#ifdef WIN32
		timeBeginPeriod(1);
		#endif
		reset();
	}

	~Clock() {
		#ifdef WIN32
		timeEndPeriod(1);
		#endif
	}

	/** Reset internal time base */
	void reset() {
		t0 = getAbs();
	}

	/**	Get system clock relative to internal time base */
	double getRel() {
		return getAbs() - t0;
	}

	/** Get absolute system clock (platform specific)
		On UNIX-like systems this will yield seconds since the epoch (1970)
		On Windows this yields seconds since the machine started (I think)
	*/
	double getAbs() {
		#ifdef WIN32
		return timeGetTime() * 0.001;
		#else
		struct timeval tv;
		gettimeofday(&tv,0);
		return (double) tv.tv_sec + (double) tv.tv_usec*1e-6;
		#endif
	}

	protected:

	double t0;
};

#endif
