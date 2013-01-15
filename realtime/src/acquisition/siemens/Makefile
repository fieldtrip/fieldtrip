# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

$(warning Compiliation appears to be broken. See bug 1532 on bugzilla.)

PLATFORM = $(shell $(CC) -dumpmachine)

# defaults, might be overwritten further down
CXXFLAGS = $(FLFLAGS) -I../../buffer/src -I../../buffer/cpp -I. -Wunused -Wall -pedantic
LDFLAGS = -lpthread ../../buffer/src/libbuffer.a
FLFLAGS = `fltk-config --cxxflags`
LDFL = `fltk-config --ldflags`
LDGL = `fltk-config --ldflags --use-gl`
EXTENSION = 
WINONLY = 

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
   FLFLAGS = -IC:/msys/1.0/local/include 
   LDFL = -mno-cygwin -lfltk -lole32 -luuid -lcomctl32 -lgdi32 -lwsock32
   LDGL = -lfltk_gl -lopengl32 -lglu32 $(LDFL)
   EXTENSION = .exe
   WINONLY = gui_streamer.exe pixeldata_to_remote_buffer.exe 
   CFLAGS  = -I../../buffer/src -I../../buffer/pthreads-win32/include -I. -Wunused -Wall -pedantic
endif

all: nii_to_buffer$(EXTENSION) opengl_client$(EXTENSION) gui_buffer_client$(EXTENSION) $(WINONLY)

###############################################################################
# Platform-independent C libraries / utilities
###############################################################################
siemensap.o: siemensap.c siemensap.h
	$(CC) -c siemensap.c -I.
	
unixtime.o: unixtime.c unixtime.h
	$(CC) -c unixtime.c $($(CXX)FLAGS)
	
###############################################################################
# Platform-independent C++ components / applications
###############################################################################
gui_buffer_client.o: gui_buffer_client.cc siemensap.h
	$(CXX) -c gui_buffer_client.cc $($(CXX)FLAGS)
	
gui_buffer_client$(EXTENSION): gui_buffer_client.o siemensap.o 
	$(CXX) -o gui_buffer_client$(EXTENSION) gui_buffer_client.o siemensap.o $(LDFLAGS) $(LDFL)

opengl_client.o: opengl_client.cc siemensap.h
	$(CXX) -c opengl_client.cc $($(CXX)FLAGS)
	
Brain3dWindow.o: Brain3dWindow.cc Brain3dWindow.h
	$(CXX) -c Brain3dWindow.cc $($(CXX)FLAGS)
	
opengl_client$(EXTENSION): opengl_client.o siemensap.o Brain3dWindow.o
	$(CXX) -o opengl_client$(EXTENSION) opengl_client.o siemensap.o Brain3dWindow.o $(LDFLAGS) $(LDGL)
	
nii_to_buffer.o: nii_to_buffer.cc
	$(CXX) -c nii_to_buffer.cc $($(CXX)FLAGS)

nii_to_buffer$(EXTENSION): nii_to_buffer.o
	$(CXX) -o nii_to_buffer$(EXTENSION) nii_to_buffer.o $(LDFLAGS)
	
save_as_nifti.o: save_as_nifti.c
	$(CC) -c save_as_nifti.c $(CFLAGS)
	
save_as_nifti$(EXTENSION): save_as_nifti.o ../../buffer/src/socketserver.o
	$(CC) -o save_as_nifti$(EXTENSION) save_as_nifti.o ../../buffer/src/socketserver.o $(LDFLAGS)

###############################################################################
# Windows specific components / acquisition software
###############################################################################
PixelDataGrabber.o: PixelDataGrabber.cc PixelDataGrabber.h
	$(CXX) -c PixelDataGrabber.cc $($(CXX)FLAGS)

gui_streamer.exe: gui_streamer.o FolderWatcher.o PixelDataGrabber.o siemensap.o
	$(CXX) -o gui_streamer gui_streamer.o FolderWatcher.o PixelDataGrabber.o siemensap.o $(LDFLAGS) $(LDFL)

pixeldata_to_remote_buffer.exe: pixeldata_to_remote_buffer.o FolderWatcher.o PixelDataGrabber.o siemensap.o
	$(CXX) -o pixeldata_to_remote_buffer pixeldata_to_remote_buffer.o FolderWatcher.o PixelDataGrabber.o siemensap.o $(LDFLAGS)

gui_streamer.o: gui_streamer.cc siemensap.h FolderWatcher.h
	$(CXX) -c gui_streamer.cc $($(CXX)FLAGS)
	
pixeldata_to_remote_buffer.o: pixeldata_to_remote_buffer.cc siemensap.h FolderWatcher.h
	$(CXX) -c pixeldata_to_remote_buffer.cc $($(CXX)FLAGS)

FolderWatcher.o: FolderWatcher.cc FolderWatcher.h
	$(CXX) -c FolderWatcher.cc $($(CXX)FLAGS)
