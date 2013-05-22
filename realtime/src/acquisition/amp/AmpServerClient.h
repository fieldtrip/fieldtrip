
#ifndef __AMP_SERVER_CLIENT_H__
#define __AMP_SERVER_CLIENT_H__

#include <string>
#include <stdint.h> 
#include <OnlineDataManager.h>
#include <FtBuffer.h>


struct AmpServerClientConfig {
	char * 	hostname;
	int		port;
	char *	amphostname;
	int		ampcommandport;
	int		ampstreamport;
	int		sfreq;
};

typedef struct {
	int64_t		ampId;
	uint64_t	length;
} AmpDataPacketHeader;

#define COMMAND_SIZE 4096
#define MAX_SAMPLES_PER_BLOCK 3000 //Guess!
#define RCV_BUFFER_SIZE 1152 * MAX_SAMPLES_PER_BLOCK
#define AMP_SAMP_SIZE 1152
#define AMP_ANALOG_CHANNELS_OFFSET 0
#define AMP_ANALOG_CHANNELS 8
#define AMP_CONSTANT_RATE 1000

#ifdef DEBUG
	#define PREFIX "DEBUG::%s:%d "
	#define DPRINTF(...) { fprintf(stdout, PREFIX ,__FILE__, __LINE__); fprintf(stdout, __VA_ARGS__); }
#else
	#define DPRINTF(...) 
#endif

#define COMMAND_FORMAT "(sendCommand cmd_%s %d %d %d)\n"

class AmpServerClient {
	public:
		AmpServerClient(struct AmpServerClientConfig * config);
		~AmpServerClient();
		
		bool connectClient();
		void disconnectClient();
		
		void start();
		void stop();
		
		void error(std::string msg);
	
		int getNumChannels();
		int getSamplingFreq();
		
		void getAmpId();
		void getAmpDetails();
		int getSubsample();
		void sendCommand(std::string cmd, int param1, int param2, int param3);
		void sendStrCommand(std::string cmd, int param1, int param2, int param3);
		
		bool getResponseInt(std::string param, unsigned int* result);
		
		int	getCurrentTime();
		
		unsigned int checkNewData();
		unsigned int readNewData(int32_t * ptr, unsigned int topass, FtEventList elist);
		
		inline void toHostOrder(char *data, int len);
		
	private:
		struct 	AmpServerClientConfig * config;
		AmpDataPacketHeader	ampDataPacketHeader;
		int 			cmdsockfd;
		int 			strsockfd;
		FILE*			strstreamin;
		FILE*			strstreamout;
		bool 			connected;
		char *	 		rcvbuffer;
		char 			r_buffer[COMMAND_SIZE];
		char			s_buffer[COMMAND_SIZE];
		unsigned int	nAmp;
		int				ampId;
		unsigned int	nchannels;
		int				lastdin1;
		int				lastdin2;
		void 			prepareCommand(std::string cmd, int param1, int param2, int param3);
		unsigned char	decodeDin(unsigned char din);
		unsigned int	subsample;
		unsigned int 	nsamples;
};




#endif
