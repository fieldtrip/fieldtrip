mex -I../src -I../../pthreads-win32/include -c buffer_gethdr.c
mex -I../src -I../../pthreads-win32/include -c buffer_getdat.c
mex -I../src -I../../pthreads-win32/include -c buffer_getevt.c
mex -I../src -I../../pthreads-win32/include -c buffer_getprp.c
mex -I../src -I../../pthreads-win32/include -c buffer_puthdr.c
mex -I../src -I../../pthreads-win32/include -c buffer_putdat.c
mex -I../src -I../../pthreads-win32/include -c buffer_putevt.c
mex -I../src -I../../pthreads-win32/include -c buffer_putprp.c
mex -I../src -I../../pthreads-win32/include -c buffer_flushhdr.c
mex -I../src -I../../pthreads-win32/include -c buffer_flushdat.c
mex -I../src -I../../pthreads-win32/include -c buffer_flushevt.c

% this is for the demo acquisition threads
mex -I../src -I../../pthreads-win32/include -c ../test/sinewave.c
mex -I../src -I../../pthreads-win32/include -c ../test/event.c

% link all the objects together into a mex file
if ispc
  % this is needed for Borland on windows, but not for other platforms
  mex -I../src -I../../pthreads-win32/include -c ../src/gettimeofday.c -I../src -I../../pthreads-win32/include
  % link all the objects together into a mex file
  mex -I../src -I../../pthreads-win32/include -L../src -L../../pthreads-win32/lib buffer.c event.obj sinewave.obj gettimeofday.obj buffer_gethdr.obj buffer_getdat.obj buffer_getevt.obj buffer_getprp.obj buffer_flushhdr.obj buffer_flushdat.obj buffer_flushevt.obj buffer_puthdr.obj buffer_putdat.obj buffer_putevt.obj buffer_putprp.obj -lbuffer -lpthreadVC2.bcb
else
  % link all the objects together into a mex file
  mex -I../src -L../src buffer.c event.o sinewave.o buffer_gethdr.o buffer_getdat.o buffer_getevt.o buffer_getprp.o buffer_flushhdr.o buffer_flushdat.o buffer_flushevt.o buffer_puthdr.o buffer_putdat.o buffer_putevt.o buffer_putprp.o -lbuffer
end

