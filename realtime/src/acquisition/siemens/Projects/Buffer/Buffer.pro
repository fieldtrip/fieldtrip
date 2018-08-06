# First declare some definitions about the directory structure
#FIELDTRIP =         ../../../../../..
#REALTIMEFOLDER=    $$FIELDTRIP/realtime
REALTIMEFOLDER =    ../../../../..
FTBUFFER =          $$REALTIMEFOLDER/src/buffer
SOURCE_DIRECTORY =  $$FTBUFFER/src
INCLUDE_DIRECTORY = $$FTBUFFER/src

CONFIG +=           static

OBJECTS_DIR =       obj

TEMPLATE =          lib

TARGET =            ../../../lib/buffer

INCLUDEPATH +=      $$FTBUFFER/src

win32:{
LIBS +=             -L$$REALTIMEFOLDER/bin/win32 -lpthreadGC2 \ #dynamic
                    -L../../lib -lwinmm -lws2_32 #static
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

