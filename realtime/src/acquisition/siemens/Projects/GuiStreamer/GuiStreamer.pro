# First declare some definitions about the directory structure
#FIELDTRIP =         ../../../../../..
#REALTIMEFOLDER=    $$FIELDTRIP/realtime
REALTIMEFOLDER=     ../../../../..
FTBUFFER =          $$REALTIMEFOLDER/src/buffer

SOURCE_DIRECTORY =  ../../src
INCLUDE_DIRECTORY = ../../include

CONFIG += static

OBJECTS_DIR =       obj

win32:{
#TARGET =            $$REALTIMEFOLDER/../bin/win32/gui_streamer
TARGET =            ../../../bin/win32/gui_streamer
}

INCLUDEPATH +=      ../../include \
                    $$FTBUFFER/src \
                    $$FTBUFFER/cpp \
                    ../../fltk #change this path to FieldTrip/external

SOURCES += \
                    $$SOURCE_DIRECTORY/gui_streamer.cc \
                    $$SOURCE_DIRECTORY\FolderWatcher.cc \
                    $$SOURCE_DIRECTORY\PixelDataGrabber.cc \
                    $$SOURCE_DIRECTORY\siemensap.c

HEADERS += \
                    $$INCLUDE_DIRECTORY\FolderWatcher.h \
                    $$INCLUDE_DIRECTORY\PixelDataGrabber.h \
                    $$INCLUDE_DIRECTORY\siemensap.h


LIBS +=             -L../../lib -lbuffer -lwinmm -lws2_32 \ #statically linked
                    -L../../lib -lfltk -lpthreadGC2 #dynamically linked, for some reason linking fltk statically didn't work
#                    -L$$REALTIMEFOLDER/bin/win32 -lpthreadGC2 #dynamically linked
