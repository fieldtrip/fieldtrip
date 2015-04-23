# First declare some definitions about the directory structure
#FIELDTRIP =         ../../../../../..
#REALTIMEFOLDER=    $$FIELDTRIP/realtime
REALTIMEFOLDER=     ../../../../..
FTBUFFER =          $$REALTIMEFOLDER/src/buffer

SOURCE_DIRECTORY =  ../../src
INCLUDE_DIRECTORY = ../../include

Debug:TARGET =      ../../../debug/nii_to_buffer
Release:TARGET =    ../../../release/nii_to_buffer

INCLUDEPATH +=      ../../include \
                    $$FTBUFFER/src \
                    $$FTBUFFER/cpp

SOURCES += \
                    $$SOURCE_DIRECTORY/nii_to_buffer.cc

LIBS +=             -L../../debug -lwinmm -lws2_32 \
                    -L../../debug -lbuffer
