# First declare some definitions about the directory structure
#FIELDTRIP =         ../../../../../..
#REALTIMEFOLDER=    $$FIELDTRIP/realtime
REALTIMEFOLDER=     ../../../../..
FTBUFFER =          $$REALTIMEFOLDER/src/buffer

SOURCE_DIRECTORY =  ../../src
INCLUDE_DIRECTORY = ../../include


Debug:TARGET =      ../../../debug/save_as_nifti
Release:TARGET =    ../../../release/save_as_nifti

INCLUDEPATH +=      ../../include \
                    $$FTBUFFER/src

SOURCES +=          $$SOURCE_DIRECTORY/save_as_nifti.c

LIBS +=             -L../../debug -lwinmm -lws2_32 \
                    -L../../debug -lbuffer
