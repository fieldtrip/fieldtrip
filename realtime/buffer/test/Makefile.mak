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
SUFFIX = _64.exe
!ELSEIF "$(PLATFORM)" == "x86"
INCPATH	 = -I../src -I../pthreads-win32/include
CFLAGS	 = $(CFLAGS) $(INCPATH) /Zi -D_X86_=1  -DWIN32 -D_WIN32 -W3 /MT
LIBS 	 = ../pthreads-win32/lib/pthreadVC2.lib ../src/libbuffer.lib ws2_32.lib mswsock.lib
LD = link
SUFFIX = .exe
!ELSE
!ERROR  Processor architecture unknown: Must specify PLATFORM=AMD64 or PLATFROM=x86
!ENDIF



all: demo

test: test_gethdr.exe test_getdat.exe test_getevt.exe test_flushhdr.exe test_flushdat.exe test_flushevt.exe test_connect.exe test_nslookup.exe test_benchmark.exe test_waitdat.exe

demo: demo_sinewave$(SUFFIX) demo_event$(SUFFIX) demo_buffer$(SUFFIX) demo_combined$(SUFFIX) demo_buffer_unix$(SUFFIX)

demo_combined$(SUFFIX): demo_combined.obj sinewave.obj 
	$(LD) /OUT:$@ $** $(LIBS) 

demo_sinewave$(SUFFIX): demo_sinewave.obj sinewave.obj
	$(LD) /OUT:$@ $** $(LIBS) 

demo_event$(SUFFIX): demo_event.obj event.obj
	$(LD) /OUT:$@ $** $(LIBS) 

demo_buffer$(SUFFIX): demo_buffer.obj 
	$(LD) /OUT:$@ $** $(LIBS) 
	
demo_buffer_unix$(SUFFIX): demo_buffer_unix.obj
	$(LD) /OUT:$@ $** $(LIBS) 

test_gethdr.exe: test_gethdr.obj
	$(LD) /OUT:$@ $** $(LIBS) 

test_getdat.exe: test_getdat.obj
	$(LD) /OUT:$@ $** $(LIBS) 

test_getevt.exe: test_getevt.obj
	$(LD) /OUT:$@ $** $(LIBS) 

test_flushhdr.exe: test_flushhdr.obj
	$(LD) /OUT:$@ $** $(LIBS) 

test_flushdat.exe: test_flushdat.obj
	$(LD) /OUT:$@ $** $(LIBS) 

test_flushevt.exe: test_flushevt.obj
	$(LD) /OUT:$@ $** $(LIBS) 

test_append.exe: test_append.obj
	$(LD) /OUT:$@ $** $(LIBS) 

test_pthread.exe: test_pthread.obj
	$(LD) /OUT:$@ $** $(LIBS) 

test_benchmark.exe: test_benchmark.obj
	$(LD) /OUT:$@ $** $(LIBS) 
	
test_nslookup.exe: test_nslookup.obj
	$(LD) /OUT:$@ $** $(LIBS) 

test_connect.exe: test_connect.obj
	$(LD) /OUT:$@ $** $(LIBS)
	
test_waitdat.exe: test_waitdat.obj
	$(LD) /OUT:$@ $** $(LIBS)
		
clean:
	del *.obj *$(SUFFIX)

distclean:
	del test_gethdr.exe test_getdat.exe test_getevt.exe test_flushhdr.exe test_flushdat.exe test_flushevt.exe
	del demo_sinewave$(SUFFIX) demo_event$(SUFFIX) demo_buffer$(SUFFIX) demo_combined$(SUFFIX)
