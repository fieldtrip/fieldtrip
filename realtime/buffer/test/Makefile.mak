##################################################
#   Visual C Makefile (NMAKE -f Makefile.mak     #
##################################################
!IF "$(PLATFORM)"==""
PLATFORM = $(PROCESSOR_ARCHITECTURE)
!ENDIF

CFLAGS = -nologo -GS
!IF "$(PLATFORM)" == "AMD64"
INCPATH	 = -I../src -I../pthreads-win64/include
CFLAGS	 = $(CFLAGS) $(INCPATH) -D_AMD64_=1 -DWIN64 -D_WIN64  -DWIN32 -D_WIN32 -W3 /MT
LIBS 	 = ../pthreads-win64/lib/pthreadVC2.lib ../src/libbuffer.lib ws2_32.lib mswsock.lib
LD = link
!ELSEIF "$(PLATFORM)" == "x86"
INCPATH	 = -I../src -I../pthreads-win32/include
CFLAGS	 = $(CFLAGS) $(INCPATH) /Zi -D_X86_=1  -DWIN32 -D_WIN32 -W3 /MT
LIBS 	 = ../pthreads-win32/lib/pthreadVC2.lib ../src/libbuffer.lib ws2_32.lib mswsock.lib
LD = link
!ELSE
!ERROR  Processor architecture unknown: Must specify PLATFORM=AMD64 or PLATFROM=x86
!ENDIF



all: demo

test: test_gethdr.exe test_getdat.exe test_getevt.exe test_flushhdr.exe test_flushdat.exe test_flushevt.exe test_connect.exe test_nslookup.exe test_benchmark.exe test_waitdat.exe

demo: demo_sinewave.exe demo_event.exe demo_buffer.exe demo_combined.exe

demo_combined.exe: demo_combined.obj sinewave.obj 
	$(LD) $** $(LIBS) 

demo_sinewave.exe: demo_sinewave.obj sinewave.obj
	$(LD) $** $(LIBS) 

demo_event.exe: demo_event.obj event.obj
	$(LD) $** $(LIBS) 

demo_buffer.exe: demo_buffer.obj 
	$(LD) $** $(LIBS) 
	
demo_buffer_unix.exe: demo_buffer_unix.obj ../src/socketserver.obj
	$(LD) $** $(LIBS) 

test_gethdr.exe: test_gethdr.obj
	$(LD) $** $(LIBS) 

test_getdat.exe: test_getdat.obj
	$(LD) $** $(LIBS) 

test_getevt.exe: test_getevt.obj
	$(LD) $** $(LIBS) 

test_flushhdr.exe: test_flushhdr.obj
	$(LD) $** $(LIBS) 

test_flushdat.exe: test_flushdat.obj
	$(LD) $** $(LIBS) 

test_flushevt.exe: test_flushevt.obj
	$(LD) $** $(LIBS) 

test_append.exe: test_append.obj
	$(LD) $** $(LIBS) 

test_pthread.exe: test_pthread.obj
	$(LD) $** $(LIBS) 

test_benchmark.exe: test_benchmark.obj
	$(LD) $** $(LIBS) 
	
test_nslookup.exe: test_nslookup.obj
	$(LD) $** $(LIBS) 

test_connect.exe: test_connect.obj
	$(LD) $** $(LIBS)
	
test_waitdat.exe: test_waitdat.obj
	$(LD) $** $(LIBS)
		
clean:
	del *.obj *.exe

distclean:
	del test_gethdr.exe test_getdat.exe test_getevt.exe test_flushhdr.exe test_flushdat.exe test_flushevt.exe
	del demo_sinewave.exe demo_event.exe demo_buffer.exe demo_combined.exe
