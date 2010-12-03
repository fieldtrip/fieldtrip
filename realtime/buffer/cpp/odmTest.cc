#include <OnlineDataManager.h>
#include <ConsoleInput.h>
#include <StringServer.h>

#define NCHAN  6
#define NBLK   237

int main() {
	ConsoleInput conIn;
	StringServer ctrlServ;
	int counter = 0;
				
	OnlineDataManager<int, float> ODM(1, NCHAN, 2000.0, GDF_INT32, DATATYPE_FLOAT32);
	int value[NCHAN];
	int speed[NCHAN];
	
	for (int i=0;i<NCHAN;i++) {
		value[i] = 0;
		speed[i] = 4*(1+i);
	}
	
	if (!ODM.useOwnServer(1972)) {
		fprintf(stderr, "Could not spawn buffer server.\n");
		return 0;
	}
	
	if (ODM.configureFromFile("config6.txt") != 0) {
		fprintf(stderr, "Configuration file is invalid\n");
		return 0;
	}
	
	ctrlServ.startListening(8000);
	
	ODM.start();
	ODM.enableSaving("test6");
   
	printf("Starting - press ESC to quit\n");
	while (1) {
		if (conIn.checkKey()) {
			int c = conIn.getKey();
			if (c==27) break; // quit
		}
		
		ctrlServ.checkRequests(ODM);
		
		int *block = ODM.provideBlock(NBLK);
		for (int j=0;j<NBLK;j++) {
			// status
			block[j*(1+NCHAN)] = 0;
			for (int i=0;i<NCHAN;i++) {
			
				value[i] += speed[i];
				if (value[i] > 1023) value[i]-=2048;
			
				block[1+i+j*(1+NCHAN)] = value[i];
			}
		}
		//printf("Block %i : %i\n", ++counter, ODM.handleBlock());
		conIn.milliSleep(100);
	}
	
	ODM.stop();
}

