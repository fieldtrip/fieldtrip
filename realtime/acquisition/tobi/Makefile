LDFLAGS = \
	-lpthread -lboost_thread -lboost_system -lboost_program_options\
	-Llib/tia/amd64 -Llib/ticpp/linux -ltia -lticpp_64 \
	-L../../buffer/src -lbuffer 

CXXFLAGS = -c -Iinclude -I ../../buffer/src

tia2ft: tia2ft.o
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: clean run

clean:
	-rm -rf *.o tia2ft
