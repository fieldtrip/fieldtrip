# The SRCS variable holds a list of all source-files for the library
SET(SRCS tcpsocket.c clientrequest.c dmarequest.c tcprequest.c util.c property.c swapbytes.c sinewave.c cleanup.c)

# If no poll-function is provided by the system, compile our own.
IF(NOT HAVE_POLL_H)
    SET(SRCS ${SRCS} win32/poll.c)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/win32) 
ENDIF(NOT HAVE_POLL_H)

# Compile and build a static library called "buffer" using the source-files in
# SRCS plus tcpserver.c
ADD_LIBRARY(buffer STATIC ${SRCS} tcpserver.c)

# Compile and build a shared library called "buffer_dll" using the same source-files
ADD_LIBRARY(buffer_dll SHARED ${SRCS} tcpserver.c)

# Shared libraries have to be linked against all libs. So we tell cmake to do this here
TARGET_LINK_LIBRARIES(buffer_dll ${PTHREAD_LIB})
TARGET_LINK_LIBRARIES(buffer_dll ${W32SOCKET_LIB})

# On Unix, we also need to link against m
IF(NOT WIN32)
    TARGET_LINK_LIBRARIES(buffer_dll m)
ENDIF(NOT WIN32)

# This sets up a "make install" target.
INSTALL(TARGETS buffer DESTINATION /usr/local/lib)
INSTALL(TARGETS buffer_dll DESTINATION /usr/local/lib)
INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/message.h DESTINATION /usr/local/include)
INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/buffer.h DESTINATION /usr/local/include)

# The sources for the standalone-buffer.
SET(STBFSRC standalone_buffer.c tcpserver.c)

# Compile and build the executable
ADD_EXECUTABLE(standalone_buffer ${STBFSRC})
# and link it against the dll. (If you want to link it against the static lib, you would have to do this:
# TARGET_LINK_LIBRARIES(standalone_buffer buffer)
# TARGET_LINK_LIBRARIES(standalone_buffer ${PTHREAD_LIB})
# TARGET_LINK_LIBRARIES(standalone_buffer ${W32SOCKET_LIB})
# IF(NOT WIN32)
#     TARGET_LINK_LIBRARIES(standalone_buffer m)
# ENDIF(NOT WIN32)
TARGET_LINK_LIBRARIES(standalone_buffer buffer_dll)

# The same for standalone_sinewave
SET(STSINSRC standalone_sinewave.c)

ADD_EXECUTABLE(standalone_sinewave ${STSINSRC})
TARGET_LINK_LIBRARIES(standalone_sinewave buffer_dll)
