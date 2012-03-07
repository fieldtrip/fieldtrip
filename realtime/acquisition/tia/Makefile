LDFLAGS = -Llib/tia/amd64 -Llib/ticpp/linux -lpthread -lboost_thread \
	-lboost_system -ltia -lticpp_64

tia2ft: tia2ft.o
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) -c -Iinclude -o $@  $^

.PHONY: clean run

run: 
	LD_LIBRARY_PATH=lib/tia/amd64 ./tia2ft

clean:
	-rm -rf *.o tia2ft
