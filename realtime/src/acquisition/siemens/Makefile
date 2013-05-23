# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
UNAME = $(shell uname)
MACHINE = $(shell uname -m)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

# defaults, might be overwritten further down
FLFLAGS = `fltk-config --cxxflags`
LDFL = `fltk-config --ldflags`
LDGL = `fltk-config --ldflags --use-gl`
CFLAGS = $(FLFLAGS) -I$(FTBUFFER)/src -I. -Wunused -Wall -pedantic
CXXFLAGS = $(CFLAGS) -I$(FTBUFFER)/cpp
LDFLAGS = -lpthread $(FTBUFFER)/src/libbuffer.a
SUFFIX = 
WINONLY = 

ifeq "$(PLATFORM)" "mingw32"
   FLFLAGS = -IC:/msys/1.0/local/include 
   LDFL = -mno-cygwin -lfltk -lole32 -luuid -lcomctl32 -lgdi32 -lwsock32
   LDGL = -lfltk_gl -lopengl32 -lglu32 $(LDFL)
   CFLAGS = -I$(FTBUFFER)/src -I$(FIELDTRIP)/realtime/src/pthreads-win32/include -I. -Wunused -Wall -pedantic
   CXXFLAGS = $(CFLAGS) -I$(FTBUFFER)/cpp
   LDFLAGS = pthreadGC2.dll $(FTBUFFER)/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
   SUFFIX = .exe
   WINONLY = gui_streamer.exe pixeldata_to_remote_buffer.exe 
endif

ifeq "$(UNAME)" "Darwin"
  # Snow Leopard 10.6 is the first Mac OS X version to ship with a 64 bit kernel, and obviously your
  # Mac has to have a 64 bit processor to utilize that kernel, but some 64 bit Mac's are not defaulting
  # to the 64 bit kernel. Furthermore, irrespective of the kernel, the gcc compiler is able to (cross)
    # compile objects, libraries and applications for 32 and 64 bit.
    ifeq "$(MACHINE)" "i386"
        # compile for 32 bit
        CFLAGS += -m32
        CXXFLAGS += -m32
        LDFLAGS += -m32
    endif
    ifeq "$(MACHINE)" "x86_64"
        # compile for 64 bit
        CFLAGS += -m64
        CXXFLAGS += -m64
        LDFLAGS += -m64
    endif
endif

all: nii_to_buffer$(SUFFIX) opengl_client$(SUFFIX) gui_buffer_client$(SUFFIX) $(WINONLY)

###############################################################################
# Platform-independent C libraries / utilities
###############################################################################
siemensap.o: siemensap.c siemensap.h
	$(CC) -c siemensap.c $(CFLAGS)

unixtime.o: unixtime.c unixtime.h
	$(CC) -c unixtime.c $(CFLAGS)

###############################################################################
# Platform-independent C++ components / applications
###############################################################################
gui_buffer_client.o: gui_buffer_client.cc siemensap.h
	$(CXX) -c gui_buffer_client.cc $(CXXFLAGS)
	
gui_buffer_client$(SUFFIX): gui_buffer_client.o siemensap.o 
	$(CXX) -o gui_buffer_client$(SUFFIX) gui_buffer_client.o siemensap.o $(LDFLAGS) $(LDFL)

opengl_client.o: opengl_client.cc siemensap.h
	$(CXX) -c opengl_client.cc $(CXXFLAGS)
	
Brain3dWindow.o: Brain3dWindow.cc Brain3dWindow.h
	$(CXX) -c Brain3dWindow.cc $(CXXFLAGS)
	
opengl_client$(SUFFIX): opengl_client.o siemensap.o Brain3dWindow.o
	$(CXX) -o opengl_client$(SUFFIX) opengl_client.o siemensap.o Brain3dWindow.o $(LDFLAGS) $(LDGL)
	
nii_to_buffer.o: nii_to_buffer.cc
	$(CXX) -c nii_to_buffer.cc $(CXXFLAGS)

nii_to_buffer$(SUFFIX): nii_to_buffer.o
	$(CXX) -o nii_to_buffer$(SUFFIX) nii_to_buffer.o $(LDFLAGS)
	
save_as_nifti.o: save_as_nifti.c
	$(CC) -c save_as_nifti.c $(CFLAGS)
	
save_as_nifti$(SUFFIX): save_as_nifti.o ../../buffer/src/socketserver.o
	$(CC) -o save_as_nifti$(SUFFIX) save_as_nifti.o ../../buffer/src/socketserver.o $(LDFLAGS)

###############################################################################
# Windows specific components / acquisition software
###############################################################################
PixelDataGrabber.o: PixelDataGrabber.cc PixelDataGrabber.h
	$(CXX) -c PixelDataGrabber.cc $(CXXFLAGS)

gui_streamer.exe: gui_streamer.o FolderWatcher.o PixelDataGrabber.o siemensap.o
	$(CXX) -o gui_streamer gui_streamer.o FolderWatcher.o PixelDataGrabber.o siemensap.o $(LDFLAGS) $(LDFL)

pixeldata_to_remote_buffer.exe: pixeldata_to_remote_buffer.o FolderWatcher.o PixelDataGrabber.o siemensap.o
	$(CXX) -o pixeldata_to_remote_buffer pixeldata_to_remote_buffer.o FolderWatcher.o PixelDataGrabber.o siemensap.o $(LDFLAGS)

gui_streamer.o: gui_streamer.cc siemensap.h FolderWatcher.h
	$(CXX) -c gui_streamer.cc $(CXXFLAGS)
	
pixeldata_to_remote_buffer.o: pixeldata_to_remote_buffer.cc siemensap.h FolderWatcher.h
	$(CXX) -c pixeldata_to_remote_buffer.cc $(CXXFLAGS)

FolderWatcher.o: FolderWatcher.cc FolderWatcher.h
	$(CXX) -c FolderWatcher.cc $(CXXFLAGS)
