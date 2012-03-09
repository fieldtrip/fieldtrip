LDFLAGS = \
	-lpthread -lboost_thread -lboost_system \
	-Llib/tia/amd64 -Llib/ticpp/linux -ltia -lticpp_64 \
	-L../../buffer/src -lbuffer

CXXFLAGS = -c -Iinclude -I ../../buffer/src

tia2ft: tia2ft.o
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@  $^

.PHONY: clean run

run: 
	LD_LIBRARY_PATH=lib/tia/amd64 ./tia2ft

clean:
	-rm -rf *.o tia2ft
