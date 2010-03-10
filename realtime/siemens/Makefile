LDFLAGS=pthreadGC2.dll ../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
CXXFLAGS=-IC:/msys/1.0/local/include -I../buffer/src -I. -g -Wall
LDFL=-lfltk -lole32 -luuid -lcomctl32 -lgdi32 -lwsock32

all: pixeldata_to_remote_buffer.exe gui_buffer_client.exe

pixeldata_to_remote_buffer.exe: pixeldata_to_remote_buffer.cc unixtime.o FolderWatcher.o
	g++ -o pixeldata_to_remote_buffer pixeldata_to_remote_buffer.cc unixtime.o FolderWatcher.o $(LDFLAGS) $(CXXFLAGS)
	
unixtime.o: unixtime.c unixtime.h
	gcc -c unixtime.c $(CXXFLAGS)
	
FolderWatcher.o: FolderWatcher.cc FolderWatcher.h
	g++ -c FolderWatcher.cc $(CXXFLAGS)
	
gui_buffer_client.exe: gui_buffer_client.cc unixtime.o
	g++ -o gui_buffer_client.exe gui_buffer_client.cc unixtime.o $(CXXFLAGS) $(LDFLAGS) $(LDFL)
