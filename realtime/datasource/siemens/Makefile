# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!
CC = gcc
CXX = g++

PLATFORM = $(shell gcc -dumpmachine)

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
   FLFLAGS = -IC:/msys/1.0/local/include 
   LDFL = -mno-cygwin -lfltk -lole32 -luuid -lcomctl32 -lgdi32 -lwsock32
   LDGL = -lfltk_gl -lopengl32 -lglu32 $(LDFL)
   SUFFIX = .exe
   WINONLY = gui_streamer.exe pixeldata_to_remote_buffer.exe 
else
   LDFLAGS = -lpthread ../../buffer/src/libbuffer.a
   FLFLAGS = `fltk-config --cxxflags`
   LDFL = `fltk-config --ldflags`
   LDGL = `fltk-config --ldflags --use-gl`
   SUFFIX = 
   WINONLY = 
endif

CXXFLAGS = $(FLFLAGS) -I../../buffer/src -I../../buffer/cpp -I. -Wunused -Wall -pedantic
 

all: nii_to_buffer$(SUFFIX) opengl_client$(SUFFIX) gui_buffer_client$(SUFFIX) $(WINONLY)

###############################################################################
# Platform-independent C libraries / utilities
###############################################################################
siemensap.o: siemensap.c siemensap.h
	gcc -c siemensap.c -I.
	
unixtime.o: unixtime.c unixtime.h
	gcc -c unixtime.c $(CXXFLAGS)
	
###############################################################################
# Platform-independent C++ components / applications
###############################################################################
gui_buffer_client.o: gui_buffer_client.cc siemensap.h
	g++ -c gui_buffer_client.cc $(CXXFLAGS)
	
gui_buffer_client$(SUFFIX): gui_buffer_client.o siemensap.o 
	g++ -o gui_buffer_client$(SUFFIX) gui_buffer_client.o siemensap.o $(LDFLAGS) $(LDFL)

opengl_client.o: opengl_client.cc siemensap.h
	g++ -c opengl_client.cc $(CXXFLAGS)
	
Brain3dWindow.o: Brain3dWindow.cc Brain3dWindow.h
	g++ -c Brain3dWindow.cc $(CXXFLAGS)
	
opengl_client$(SUFFIX): opengl_client.o siemensap.o Brain3dWindow.o
	g++ -o opengl_client$(SUFFIX) opengl_client.o siemensap.o Brain3dWindow.o $(LDFLAGS) $(LDGL)
	
nii_to_buffer.o: nii_to_buffer.cc
	g++ -c nii_to_buffer.cc $(CXXFLAGS)

nii_to_buffer$(SUFFIX): nii_to_buffer.o
	g++ -o nii_to_buffer$(SUFFIX) nii_to_buffer.o $(LDFLAGS)

###############################################################################
# Windows specific components / acquisition software
###############################################################################
PixelDataGrabber.o: PixelDataGrabber.cc PixelDataGrabber.h
	g++ -c PixelDataGrabber.cc $(CXXFLAGS)

gui_streamer.exe: gui_streamer.o FolderWatcher.o PixelDataGrabber.o siemensap.o
	g++ -o gui_streamer gui_streamer.o FolderWatcher.o PixelDataGrabber.o siemensap.o $(LDFLAGS) $(LDFL)

pixeldata_to_remote_buffer.exe: pixeldata_to_remote_buffer.o FolderWatcher.o PixelDataGrabber.o siemensap.o
	g++ -o pixeldata_to_remote_buffer pixeldata_to_remote_buffer.o FolderWatcher.o PixelDataGrabber.o siemensap.o $(LDFLAGS)

gui_streamer.o: gui_streamer.cc siemensap.h FolderWatcher.h
	g++ -c gui_streamer.cc $(CXXFLAGS)
	
pixeldata_to_remote_buffer.o: pixeldata_to_remote_buffer.cc siemensap.h FolderWatcher.h
	g++ -c pixeldata_to_remote_buffer.cc $(CXXFLAGS)

FolderWatcher.o: FolderWatcher.cc FolderWatcher.h
	g++ -c FolderWatcher.cc $(CXXFLAGS)
