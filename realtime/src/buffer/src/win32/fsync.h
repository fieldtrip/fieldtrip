#ifndef _FSYNC
#define _FSYNC

/*
 * work around lack of fsync on MinGW, see
 * https://www.mail-archive.com/bug-gnulib@gnu.org/msg11472.html
*/

int fsync(int fd);

#endif
