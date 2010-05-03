##################################################
#   Visual C Makefile (NMAKE -f Makefile.mak)    #
##################################################
!IF "$(PLATFORM)"==""
PLATFORM = $(PROCESSOR_ARCHITECTURE)
!ENDIF

CFLAGS = -nologo -GS
!IF "$(PLATFORM)" == "AMD64"
INCPATH	 = -I. -I../pthreads-win64/include
CFLAGS	 = $(CFLAGS) $(INCPATH) /Zi -D_AMD64_=1 -DWIN64 -D_WIN64  -DWIN32 -D_WIN32 -W3
LIBFLAGS = /MACHINE:X64
!ELSEIF "$(PLATFORM)" == "x64"
# USE THIS FOR CROSS-COMPILING!
INCPATH	 = -I. -I../pthreads-win64/include
CFLAGS	 = $(CFLAGS) $(INCPATH) /Zi -D_AMD64_=1 -DWIN64 -D_WIN64  -DWIN32 -D_WIN32 -W3
LIBFLAGS = /MACHINE:X64 /NODEFAULTLIB:libcmt
CC = "C:\Program Files\Microsoft Visual Studio 9.0\VC\bin\x86_amd64\cl.exe"
!ELSEIF "$(PLATFORM)" == "x86"
INCPATH	 = -I. -I../pthreads-win32/include
CFLAGS	 = $(CFLAGS) $(INCPATH) /Zi -D_X86_=1  -DWIN32 -D_WIN32 -W3
LIBFLAGS = /MACHINE:X86
!ELSE
!ERROR  Processor architecture unknown: Must specify PLATFORM=AMD64 or PLATFROM=x86 or PLATFORM=x64
!ENDIF

all: libbuffer.lib

libbuffer.lib: tcpserver.obj tcpsocket.obj tcprequest.obj clientrequest.obj dmarequest.obj cleanup.obj util.obj printstruct.obj swapbytes.obj property.obj gettimeofday.obj
	lib $(LIBFLAGS) /OUT:libbuffer.lib $**

gettimeofday.obj: win32/gettimeofday.c
	$(CC) /c $(CFLAGS) win32/gettimeofday.c
	
%.obj: %.c buffer.h message.h swapbytes.h socket_includes.h unix_includes.h
	$(CC) /c $(CFLAGS) $*.c

clean:
	del *.obj *.lib
