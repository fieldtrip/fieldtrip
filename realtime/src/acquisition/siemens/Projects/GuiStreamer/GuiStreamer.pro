# First declare some definitions about the directory structure
#FIELDTRIP =         ../../../../../..
#REALTIMEFOLDER=    $$FIELDTRIP/realtime
REALTIMEFOLDER=     ../../../../..
FTBUFFER =          $$REALTIMEFOLDER/src/buffer

CONFIG += static

SOURCE_DIRECTORY =  ../../src
INCLUDE_DIRECTORY = ../../include

Debug:TARGET =      ../../../debug/gui_streamer
Release:TARGET =    ../../../release/gui_streamer

INCLUDEPATH +=      ../../include \
                    $$FTBUFFER/src \
                    $$FTBUFFER/cpp \
                    ../../fltk

SOURCES += \
                    $$SOURCE_DIRECTORY/gui_streamer.cc \
                    $$SOURCE_DIRECTORY\FolderWatcher.cc \
                    $$SOURCE_DIRECTORY\PixelDataGrabber.cc \
                    $$SOURCE_DIRECTORY\siemensap.c

HEADERS += \
                    $$INCLUDE_DIRECTORY\FolderWatcher.h \
                    $$INCLUDE_DIRECTORY\PixelDataGrabber.h \
                    $$INCLUDE_DIRECTORY\siemensap.h


Debug:LIBS +=       -L../../debug -lwinmm -lws2_32 -lbuffer -lfltk
Release:LIBS +=     -L../../release -lwinmm -lws2_32 -lbuffer -lfltk
