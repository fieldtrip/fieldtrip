LDFLAGS = \
	-lpthread -lboost_thread -lboost_system -lboost_program_options\
	-Llib/tia/amd64 -Llib/ticpp/linux -ltia -lticpp_64 \
	-L../../buffer/src -lbuffer 

CXXFLAGS_DEBUG = -g -O0
CXXFLAGS = -c -Iinclude -I ../../buffer/src $(CXXFLAGS_DEBUG)

tia2ft: tia2ft.o
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: clean run

clean:
	-rm -rf *.o tia2ft
