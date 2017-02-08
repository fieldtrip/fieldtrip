/** Simple class wrapping up functions for checking key presses in a console,
    tested on Windows, OS X and Linux.

    (C) 2010 Stefan Klanke
 */

#ifndef __ConsoleInput_h
#define __ConsoleInput_h

#ifdef WIN32
	#include <conio.h>
	#include <windows.h>
#else
	#include <unistd.h>
	#include <sys/ioctl.h>
	#include <termios.h>
#endif

class ConsoleInput {
	public:

	ConsoleInput() {
		#ifdef WIN32
		#else
			termios term;

			tcgetattr(0, &oldTerm);
			term = oldTerm;
			term.c_lflag &= ~ICANON;
			tcsetattr(0, TCSANOW, &term);
			setbuf(stdin, NULL);
		#endif
    }

	~ConsoleInput() {
		#ifdef WIN32
		#else
			tcsetattr(0, TCSANOW, &oldTerm);
		#endif
	};

	bool checkKey() {
		#ifdef WIN32
			return _kbhit();
		#else
			int bytesWaiting;
			ioctl(0, FIONREAD, &bytesWaiting);
			return bytesWaiting > 0;
		#endif
	}

	int getKey() {
		#ifdef WIN32
			return getch();
		#else
			return getchar();
		#endif
	}

	void milliSleep(int ms) {
		#ifdef WIN32
			Sleep(ms);
		#else
			usleep(ms*1000);
		#endif
	}

	protected:
	#ifdef WIN32
	#else
	termios oldTerm;
	#endif
};

#endif
