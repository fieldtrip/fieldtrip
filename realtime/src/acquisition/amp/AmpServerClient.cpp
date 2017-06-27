/* Copyright (C) 2013 Federico Raimondo
 * Applied Artificial Intelligence Lab
 * Computer Sciences Department
 * University of Buenos Aires, Argentina
 *
 * This file is part of Amp2ft
 *
 * Amp2ft is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Amp2ft is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Amp2ft.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>

#include "platform_includes.h"
#ifndef COMPILER_MINGW
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#endif

#include <AmpServerClient.h>


AmpServerClient::AmpServerClient(struct AmpServerClientConfig * config) {
	this->config = config;
	this->connected = false;
	this->rcvbuffer = NULL;

	/* Get 1 sample each every 1000/sfreq samples
	 * cause AmpServer always sends 1000 samples per sec
	 */
	 this->subsample = AMP_CONSTANT_RATE / config->sfreq;
	 this->nsamples=0;
}

AmpServerClient::~AmpServerClient() {
	if (this->connected) {
		this->disconnectClient();
	}
	if (this->rcvbuffer != NULL) {
		free((void*)this->rcvbuffer);
	}
}

bool AmpServerClient::connectClient() {
	struct sockaddr_in serv_addr;
	struct hostent *server;
	this->cmdsockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (this->cmdsockfd < 0)
        error("ERROR opening socket");
    server = gethostbyname(this->config->amphostname);
    if (server == NULL) {
        error("ERROR, no such host\n");
    }
    bzero((char *) &serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    bcopy((char *)server->h_addr,
         (char *)&serv_addr.sin_addr.s_addr,
         server->h_length);
    serv_addr.sin_port = htons(this->config->ampcommandport);
    if (connect(this->cmdsockfd,(struct sockaddr *) &serv_addr,sizeof(serv_addr)) < 0)
        error("ERROR connecting");


	this->strsockfd = socket(AF_INET, SOCK_STREAM, 0);
	serv_addr.sin_port = htons(this->config->ampstreamport);
    if (connect(this->strsockfd,(struct sockaddr *) &serv_addr,sizeof(serv_addr)) < 0)
        error("ERROR connecting");

    this->strstreamin = fdopen(this->strsockfd, "r");
    this->strstreamout = fdopen(this->strsockfd, "w");



	/* Initialize commands */
	sendCommand("SetPower", 0, 0, 1);
	usleep(3000);


	getAmpId();

	getAmpDetails();

	this->rcvbuffer = (char*)malloc(RCV_BUFFER_SIZE);
	if (this->rcvbuffer == NULL) error("ERROR not enough memory for receive buffer");
	DPRINTF("Receive buffer at %p\n", this->rcvbuffer);
	memset(r_buffer, 0, COMMAND_SIZE);


	sendCommand("SetDecimatedRate", 0, 0, this->config->sfreq);

	return true;
}

int AmpServerClient::getSubsample() {
	return this->subsample;
}

void AmpServerClient::disconnectClient() {
	close(this->strsockfd);
	close(this->cmdsockfd);
	this->connected = false;

}
void AmpServerClient::start() {
	sendCommand("SetPower", this->ampId, 0, 1);
	sendCommand("Start", this->ampId, 0, 0);
	sendCommand("DefaultAcquisitionState", 0, 0, 0);
	sendStrCommand("ListenToAmp", this->ampId, 0, 0);
	usleep(3000);
}

void AmpServerClient::stop() {
	sendStrCommand("StopListeningToAmp", this->ampId, 0, 0);
	sendCommand("Stop", this->ampId, 0, 0);
	sendCommand("DefaultAcquisitionState", 0, 0, 0);
	sendCommand("SetPower", this->ampId, 0, 0);
}

void AmpServerClient::error(std::string msg) {
    perror(msg.c_str());
    exit(0);
}

void AmpServerClient::sendCommand(std::string cmd, int param1, int param2, int param3) {
	DPRINTF("Sending command %s\n", cmd.c_str());
	int n;
	this->prepareCommand(cmd, param1, param2, param3);
	DPRINTF("Sending command %s => %s", cmd.c_str(), this->s_buffer);
	n = write(this->cmdsockfd, this->s_buffer, strlen(this->s_buffer));
	if (n < 0)
		error("ERROR writing to socket");

	fsync(this->strsockfd);

	/*do {
		n = read(this->cmdsockfd, &r_buffer[l], 1);
		if (n > 0) {
			l++;
		}
	} while (n == 0 || n > 0 && r_buffer[l-1] != '\0' && r_buffer[l-1] != '\n');
	r_buffer[l>1?l-1:l]='\0';*/
	n = read(this->cmdsockfd, r_buffer, COMMAND_SIZE);
	if (n < 0) {
		error("ERROR reading response");
	}
	r_buffer[n>1?n-1:0] = '\0';
	DPRINTF("Received %d bytes in buffer\nContent >>>\n%s \n<<< End Content\n", n, r_buffer);
}


void AmpServerClient::sendStrCommand(std::string cmd, int param1, int param2, int param3) {
	DPRINTF("Sending stream command %s\n", cmd.c_str());

	int n;
	this->prepareCommand(cmd, param1, param2, param3);
	DPRINTF("Sending command %s => %s", cmd.c_str(), this->s_buffer);
	n = fwrite(this->s_buffer, sizeof(char), strlen(this->s_buffer), this->strstreamout);
	if (n < 0)
		error("ERROR writing to socket");

	fflush(this->strstreamout);

}

void AmpServerClient::getAmpId() {
	DPRINTF("Getting number of amps\n");
	sendCommand("NumberOfAmps", 0, 0, 0);

	if (! getResponseInt("number_of_amps", &this->nAmp)) {
		error("ERROR Cannot get number of amps");
	}

	DPRINTF("Number of amps detected = %d\n", this->nAmp);
	this->ampId = this->nAmp-1;

}

void AmpServerClient::getAmpDetails() {
	DPRINTF("Getting details of amp\n");
	sendCommand("GetAmpDetails", this->ampId, 0, 0);
	if (! getResponseInt("number_of_channels", &this->nchannels)) {
		error("ERROR Cannot get number of channels");
	}
	this->nchannels += AMP_ANALOG_CHANNELS;
}

int AmpServerClient::getNumChannels() {
	return this->nchannels;
}

int AmpServerClient::getSamplingFreq() {
	return config->sfreq;
}
unsigned int AmpServerClient::checkNewData() {
	//DPRINTF("Checking for new data\n");
	if (fread((char *)&(this->ampDataPacketHeader), sizeof(AmpDataPacketHeader), 1, this->strstreamin) == 1) {
		this->toHostOrder((char*)&(this->ampDataPacketHeader), 2);
		//DPRINTF("New data from amp %lu (%lu bytes %lu samples)\n", \
			this->ampDataPacketHeader.ampId, this->ampDataPacketHeader.length, this->ampDataPacketHeader.length/AMP_SAMP_SIZE);
		return this->ampDataPacketHeader.length/AMP_SAMP_SIZE;
	}
	return 0;
}

unsigned int AmpServerClient::readNewData(int32_t * ptr, unsigned int topass, FtEventList elist) {
	//DPRINTF("Reading new data\n");
	uint64_t size = this->ampDataPacketHeader.length;
	size_t readed;
	unsigned int retorno = 0;
	unsigned int vsamp = 0;
	unsigned int rsamp = 0;
	unsigned int pcount = topass;
	unsigned int startsamp = this->nsamples;
	while (startsamp % this->subsample != 0) {
		startsamp++;
		rsamp++;
	}
	unsigned char din1;
	unsigned char din2;
	if (size > 0) {
		if ((readed = fread(this->rcvbuffer, size, 1, this->strstreamin)) > 0) {
			//~ DPRINTF("Got %lu bytes (%lu samples) starting at sample %d\n", size, size/AMP_SAMP_SIZE, rsamp);
			while (rsamp < size/AMP_SAMP_SIZE && pcount > 0) {
				//DPRINTF("Processing real sample %d (v sample %d)\n", rsamp, vsamp);
				int offset = (rsamp * AMP_SAMP_SIZE);
				din1 = decodeDin(this->rcvbuffer[offset+24]);
				din2 = decodeDin(this->rcvbuffer[offset+25]);
				if (this->lastdin1 != din1) {
					elist.add(vsamp, "DIN1", din1);
					this->lastdin1 = din1;
				}
				if (this->lastdin2 != din2) {
					elist.add(vsamp, "DIN2", din2);
					this->lastdin2 = din2;
				}
				offset+= 32;
				for (unsigned int chan = 0; chan < this->nchannels; chan++) {
					//~ DPRINTF("Processing channel %d\n", chan);
					int32_t * p = reinterpret_cast<int32_t*>( this->rcvbuffer + offset);
					*p = ntohl( *p );
					ptr[vsamp * this->nchannels + chan] = *p;
					offset += sizeof(int);

				}
				vsamp++;
				rsamp += this->subsample;
				pcount --;
			}
			this->nsamples += size/AMP_SAMP_SIZE;
			retorno = vsamp;
		} else {
			DPRINTF("Nothing to read! (%lu of %lu)\n", readed, size)
		}

	}

	return retorno;
}


inline void AmpServerClient::toHostOrder(char *data, int len) {
	int32_t *tmpData = reinterpret_cast<int32_t *>(data);
	int location;
	int32_t tmpValue;
	for(int i=0;i<len; i++) {
		location = i*2;
		tmpValue = -1;
		tmpData[location] =ntohl(tmpData[location]);
		tmpData[location+1] =ntohl(tmpData[location+1]);
		tmpValue = tmpData[location];
		tmpData[location] = tmpData[location+1];
		tmpData[location+1] = tmpValue;
	}
}



int AmpServerClient::getCurrentTime() {
	unsigned int result;
	DPRINTF("Getting current time\n");
	sendCommand("GetCurrentTime", this->ampId, 0, 0);
	if (! getResponseInt("current_time", &result)) {
		error("ERROR Cannot get current time");
	}
	return result;
}

void AmpServerClient::prepareCommand(std::string cmd, int param1, int param2, int param3) {
	sprintf(this->s_buffer, COMMAND_FORMAT, cmd.c_str(), param1, param2, param3);
}

bool AmpServerClient::getResponseInt(std::string param, unsigned int* result) {
	const char * mystring = param.c_str();
	char * start = strstr(r_buffer, mystring);
	if ( start != NULL) {
		start += strlen(mystring) +1;
		size_t end = strcspn(start, ")");
		start[end] = 0;
		*result = atoi(start);
		return true;
	}
	return false;
}

unsigned char AmpServerClient::decodeDin(unsigned char din) {
	unsigned char result = 0;
	if ((din & 1) != 0) {
		result += 64;
	}
	if ((din & 2) != 0) {
		result += 2;
	}
	if ((din & 4) != 0) {
		result += 1;
	}
	if ((din & 8) != 0) {
		result += 0;
	}
	if ((din & 16) != 0) {
		result += 16;
	}
	if ((din & 32) != 0) {
		result += 4;
	}
	if ((din & 64) != 0) {
		result += 8;
	}
	if ((din & 128) != 0) {
		result += 0;
	}
	return result;


}
