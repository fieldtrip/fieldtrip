LIBS = -Llib/tia/amd64 -Llib/ticpp/linux -lpthread -lboost_thread \
	-lboost_system -ltia -lticpp_64
LDFLAGS = $(LIBS)

run: tia_client_main
	LD_LIBRARY_PATH=lib/tia/amd64 ./tia_client_main

tia_client_main: tia_client_main.o
	$(CXX) $(LDFLAGS)  -o $@ $^

%.o: %.cpp
	$(CXX) -c -Iinclude -o $@  $^

.PHONY: clean run
clean:
	-rm -rf *.o tia_client_main
