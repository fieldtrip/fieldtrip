LDFLAGS=pthreadGC2.dll ../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
CXXFLAGS=-IC:/msys/1.0/local/include -I../buffer/src -I. -g -Wall
LDFL=-lfltk -lole32 -luuid -lcomctl32 -lgdi32 -lwsock32

all: pixeldata_to_remote_buffer.exe gui_buffer_client.exe gui_streamer.exe

gui_streamer.exe: gui_streamer.o FolderWatcher.o PixelDataGrabber.o siemensap.o
	g++ -o gui_streamer gui_streamer.o FolderWatcher.o PixelDataGrabber.o siemensap.o $(LDFLAGS) $(LDFL)

pixeldata_to_remote_buffer.exe: pixeldata_to_remote_buffer.o FolderWatcher.o PixelDataGrabber.o siemensap.o
	g++ -o pixeldata_to_remote_buffer pixeldata_to_remote_buffer.o FolderWatcher.o PixelDataGrabber.o siemensap.o $(LDFLAGS)
	
pixeldata_to_remote_buffer.o: pixeldata_to_remote_buffer.cc
	g++ -c pixeldata_to_remote_buffer.cc $(CXXFLAGS)

gui_streamer.o: gui_streamer.cc
	g++ -c gui_streamer.cc $(CXXFLAGS)
	
unixtime.o: unixtime.c unixtime.h
	gcc -c unixtime.c $(CXXFLAGS)
	
FolderWatcher.o: FolderWatcher.cc FolderWatcher.h
	g++ -c FolderWatcher.cc $(CXXFLAGS)
	
gui_buffer_client.o: gui_buffer_client.cc 	
	g++ -c gui_buffer_client.cc $(CXXFLAGS)
	
gui_buffer_client.exe: gui_buffer_client.o siemensap.o 
	g++ -o gui_buffer_client.exe gui_buffer_client.o siemensap.o $(LDFLAGS) $(LDFL)

siemensap.o: siemensap.c siemensap.h
	gcc -c siemensap.c -I.

PixelDataGrabber.o: PixelDataGrabber.cc PixelDataGrabber.h
	g++ -c PixelDataGrabber.cc $(CXXFLAGS)