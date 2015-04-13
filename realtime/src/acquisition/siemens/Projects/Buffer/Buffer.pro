# First declare some definitions about the directory structure
#FIELDTRIP =         ../../../../../..
#REALTIMEFOLDER=    $$FIELDTRIP/realtime
REALTIMEFOLDER =    ../../../../..
FTBUFFER =          $$REALTIMEFOLDER/src/buffer
SOURCE_DIRECTORY =  $$FTBUFFER/src
INCLUDE_DIRECTORY = $$FTBUFFER/src

TEMPLATE =          lib

Debug:TARGET =      ../../../debug/buffer
Release:TARGET =    ../../../release/buffer

INCLUDEPATH +=      $$FTBUFFER/src

win32:{
Debug:LIBS +=       -L../../debug -lpthreadGC2 \
                    -L../../debug -lws2_32

Release:LIBS +=     -L../../release -lpthreadGC2 \
                    -L../../release -lws2_32
}

HEADERS += \
                    $$SOURCE_DIRECTORY/buffer.h \
                    $$SOURCE_DIRECTORY/extern.h \
                    $$SOURCE_DIRECTORY/message.h \
                    $$SOURCE_DIRECTORY/platform.h \
                    $$SOURCE_DIRECTORY/platform_includes.h \
                    $$SOURCE_DIRECTORY/rdadefs.h \
                    $$SOURCE_DIRECTORY/rdaserver.h \
                    $$SOURCE_DIRECTORY/socketserver.h \
                    $$SOURCE_DIRECTORY/swapbytes.h \
                    $$SOURCE_DIRECTORY/tiaserver.h

SOURCES += \
                    $$INCLUDE_DIRECTORY/cleanup.c \
                    $$INCLUDE_DIRECTORY/clientrequest.c \
                    $$INCLUDE_DIRECTORY/dmarequest.c \
                    $$INCLUDE_DIRECTORY/endianutil.c \
                    $$INCLUDE_DIRECTORY/extern.c \
                    $$INCLUDE_DIRECTORY/printstruct.c \
                    $$INCLUDE_DIRECTORY/rdaserver.c \
                    $$INCLUDE_DIRECTORY/socketserver.c \
                    $$INCLUDE_DIRECTORY/swapbytes.c \
                    $$INCLUDE_DIRECTORY/tcprequest.c \
                    $$INCLUDE_DIRECTORY/tcpserver.c \
                    $$INCLUDE_DIRECTORY/tcpsocket.c \
                    $$INCLUDE_DIRECTORY/tiaserver.c \
                    $$INCLUDE_DIRECTORY/util.c

