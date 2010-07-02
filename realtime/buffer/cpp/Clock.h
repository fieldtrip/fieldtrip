#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

class Clock {
	public:
	
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
	
	void reset() {
		t0 = getAbs();
	}
	
	double getRel() {
		return getAbs() - t0;
	}
	
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

