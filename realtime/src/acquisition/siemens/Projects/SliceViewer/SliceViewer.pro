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
                    $$SOURCE_DIRECTORY/opengl_client.cc
                    $$SOURCE_DIRECTORY/Brain3dWindow.cc

HEADERS += \
                    $$INCLUDE_DIRECTORY/Brain3dWindow.cc

Debug:LIBS +=       -L../../debug -lwinmm -lws2_32 -lbuffer -lfltk
Release:LIBS +=     -L../../release -lwinmm -lws2_32 -lbuffer -lfltk
