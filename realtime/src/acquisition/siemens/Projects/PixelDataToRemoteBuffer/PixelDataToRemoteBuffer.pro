# First declare some definitions about the directory structure
FIELDTRIP =     ../../../../../..
#REALTIMEFOLDER=realtime
# First declare some definitions about the directory structure
#FIELDTRIP =         ../../../../../..
#REALTIMEFOLDER=    $$FIELDTRIP/realtime
REALTIMEFOLDER=     ../../../../..
FTBUFFER =          $$REALTIMEFOLDER/src/buffer

SOURCE_DIRECTORY =  ../../src
INCLUDE_DIRECTORY = ../../include

Debug:TARGET =      ../../../debug/pixeldata_to_remote_buffer
Release:TARGET =    ../../../release/pixeldata_to_remote_buffer

INCLUDEPATH +=      ../../include \
                    $$FTBUFFER/src \
                    $$FTBUFFER/cpp

HEADERS += \
                    $$INCLUDE_DIRECTORY/FolderWatcher.h \
                    $$INCLUDE_DIRECTORY/PixelDataGrabber.h \
                    $$INCLUDE_DIRECTORY/siemensap.h

SOURCES += \
                    $$SOURCE_DIRECTORY/FolderWatcher.cc \
                    $$SOURCE_DIRECTORY/PixelDataGrabber.cc \
                    $$SOURCE_DIRECTORY/siemensap.c \
                    $$SOURCE_DIRECTORY/pixeldata_to_remote_buffer.cc

LIBS +=             -L../../debug -lwinmm -lws2_32 \
                    -L../../debug -lbuffer


