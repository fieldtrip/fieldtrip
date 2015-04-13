# First declare some definitions about the directory structure
#FIELDTRIP =         ../../../../../..
#REALTIMEFOLDER=    $$FIELDTRIP/realtime
REALTIMEFOLDER=     ../../../../..
FTBUFFER =          $$REALTIMEFOLDER/src/buffer

SOURCE_DIRECTORY =  ../../src
INCLUDE_DIRECTORY = ../../include

Debug:TARGET =      ../../../debug/gui_buffer_client
Release:TARGET =    ../../../release/gui_buffer_client

INCLUDEPATH +=      ../../include \
                    $$FTBUFFER/src \
                    $$FTBUFFER/cpp \
                    ../../fltk

SOURCES += \
                    $$SOURCE_DIRECTORY/gui_buffer_client.cc \
                    $$SOURCE_DIRECTORY\siemensap.c

HEADERS += \
                    $$INCLUDE_DIRECTORY\siemensap.h

LIBS +=             -L../../debug -lbuffer \
                    -L../../fltk/lib -lfltk

