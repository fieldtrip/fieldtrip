# Set the projects name (Must be set, can be changed to whatever...)
PROJECT(FCDCBuffer)

# Include Macro to check whether header files can be included
INCLUDE(CheckIncludeFileCXX)

# Add compiler-flag to produce debug-symbols
ADD_DEFINITIONS(-g)

# Now we look whether Pthreads are installed. Note that the name of the lib depends on
# the system we compile for.
MESSAGE(STATUS "Looking for pthread library...")
IF(WIN32)
    FIND_LIBRARY(PTHREAD_LIB NAMES pthreadGC2 pthread PATHS /usr/lib /usr/local/lib)
ELSE(WIN32)
    FIND_LIBRARY(PTHREAD_LIB NAMES pthread PATHS /usr/lib /usr/local/lib)
ENDIF(WIN32)

IF(PTHREAD_LIB)
    MESSAGE(STATUS "Found pthread library: " ${PTHREAD_LIB})
ELSE(PTHREAD_LIB)
    MESSAGE(FATAL_ERROR "Pthread library not found...")
ENDIF(PTHREAD_LIB)

# Check for existense and usability of Pthreads headers
CHECK_INCLUDE_FILE_CXX(pthread.h HAVE_PTHREAD_H)
IF(NOT HAVE_PTHREAD_H)
    MESSAGE(FATAL_ERROR "Please check your Pthread setup!")
ENDIF(NOT HAVE_PTHREAD_H)

CHECK_INCLUDE_FILE_CXX(sched.h HAVE_SCHED_H)
IF(NOT HAVE_SCHED_H)
    MESSAGE(FATAL_ERROR "Please check your Pthread setup!")
ENDIF(NOT HAVE_SCHED_H)

# MinGW does not provide a poll.h. So if it is not found, we will use one that is
# provided by us.
CHECK_INCLUDE_FILE_CXX(poll.h HAVE_POLL_H)
IF(NOT HAVE_POLL_H)
    MESSAGE(STATUS "Using poll.h emulation")
ENDIF(NOT HAVE_POLL_H)

# If we compile for windows, we have to link against the winsock api...
IF(WIN32)
    MESSAGE(STATUS "Looking for win32 socket api...")
    #FIND_LIBRARY(W32SOCKET_LIB wsock32)
    FIND_LIBRARY(W32SOCKET_LIB ws2_32)
    IF(NOT W32SOCKET_LIB)
        MESSAGE(FATAL_ERROR "Not found")
    ENDIF(NOT W32SOCKET_LIB)
ELSE(WIN32)
    SET(W32SOCKET_LIB "")
ENDIF(WIN32)


# Decent into the src folder (and process the CMakeLists.txt there)
SUBDIRS(src)


