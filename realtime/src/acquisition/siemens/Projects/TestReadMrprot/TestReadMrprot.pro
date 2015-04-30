SOURCE_DIRECTORY =  ../../src
INCLUDE_DIRECTORY = ../../include

Debug:TARGET =      test
Release:TARGET =    test

INCLUDEPATH +=      ../../include

HEADERS += \
                    $$INCLUDE_DIRECTORY/siemensap.h

SOURCES +=          $$SOURCE_DIRECTORY/tsap.c \
                    $$SOURCE_DIRECTORY/siemensap.c
