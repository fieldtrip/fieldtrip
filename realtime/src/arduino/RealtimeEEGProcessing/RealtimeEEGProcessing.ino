
// This is part of a realtime EEG processing demonstration. 
//
// This particular device connnects to the FieldTrip buffer over Ethernet, 
// reads a short snippet of EEG data, computes the RMS value as a measure 
// of overall signal intensity and passes this RMS value on to a wireless 
// RFM12b (433/886 MHz) connection for remote feedback.
//
// It is based on a Sparkfun Pro Micro Arduino Pro Mini 3V3. 
// Ethernet is provided by a WIZnet W5100 network module. 
// Wireless connectivity is provided by a RFM12B module.
//
// See http://www.fieldtriptoolbox.org/development/realtime/arduino for documentation.

// Copyright (C) 2013, Robert Oostenveld, DCCN
//
// This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
// for the documentation and details.
//
//    FieldTrip is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    FieldTrip is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
//
// $Id$

#include <SPI.h>
#include <Ethernet.h>
#include <RF12.h>
#include <JeeLib.h>

/***************************************************************************************************
 * FieldTrip buffer 
 **************************************************************************************************/

// define the version of the message packet 
#define VERSION    (uint16_t)0x0001

// these define the commands that can be used, which are split over the two available bytes 
#define PUT_HDR    (uint16_t)0x0101
#define PUT_DAT    (uint16_t)0x0102
#define PUT_EVT    (uint16_t)0x0103
#define PUT_OK     (uint16_t)0x0104
#define PUT_ERR    (uint16_t)0x0105

#define GET_HDR    (uint16_t)0x0201 
#define GET_DAT    (uint16_t)0x0202
#define GET_EVT    (uint16_t)0x0203
#define GET_OK     (uint16_t)0x0204
#define GET_ERR    (uint16_t)0x0205

#define FLUSH_HDR  (uint16_t)0x0301 
#define FLUSH_DAT  (uint16_t)0x0302
#define FLUSH_EVT  (uint16_t)0x0303
#define FLUSH_OK   (uint16_t)0x0304
#define FLUSH_ERR  (uint16_t)0x0305

#define WAIT_DAT   (uint16_t)0x0402
#define WAIT_OK    (uint16_t)0x0404
#define WAIT_ERR   (uint16_t)0x0405

// these are used in the data_t and event_t structure 
#define DATATYPE_CHAR    (uint32_t)0
#define DATATYPE_UINT8   (uint32_t)1
#define DATATYPE_UINT16  (uint32_t)2
#define DATATYPE_UINT32  (uint32_t)3
#define DATATYPE_UINT64  (uint32_t)4
#define DATATYPE_INT8    (uint32_t)5
#define DATATYPE_INT16   (uint32_t)6
#define DATATYPE_INT32   (uint32_t)7
#define DATATYPE_INT64   (uint32_t)8
#define DATATYPE_FLOAT32 (uint32_t)9
#define DATATYPE_FLOAT64 (uint32_t)10

// a packet that is sent over the network (or to disk) should contain the following 
typedef struct {
  uint16_t version;   // see VERSION 
  uint16_t command;   // see PUT_xxx, GET_xxx and FLUSH_xxx 
  uint32_t bufsize;   // size of the buffer in bytes 
} 
messagedef_t; // 8 bytes

// the header definition is fixed, except for the channel labels and other chunks 
typedef struct {
  uint32_t  nchans;
  uint32_t  nsamples;
  uint32_t  nevents;
  float     fsample;
  uint32_t  data_type;
  uint32_t  bufsize;   // size of the buffer in bytes 
} 
headerdef_t; // 24 bytes

// the data definition is fixed
typedef struct {
  uint32_t nchans;
  uint32_t nsamples;
  uint32_t data_type;
  uint32_t bufsize;   // size of the buffer in bytes 
} 
datadef_t; // 16 bytes

// prepare the GET_HDR and GET_DAT requests
byte get_hdr[8] = {
  0x00, 0x01, 0x02, 0x01, 0x00, 0x00, 0x00, 0x00};
byte get_dat[16] = {
  0x00, 0x01, 0x02, 0x02, 0x00, 0x00, 0x00, 0x08, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}; // the last 2*4 bytes will be updated further down
byte buf[256];

messagedef_t request, response;
headerdef_t  header;
datadef_t    data;

/***************************************************************************************************
 * datatype conversion 
 **************************************************************************************************/

union {
  byte  bval[4];
  float fval;
} 
single; // for conversion between bytes and float

// the following is a collection of macro's to deal with converting the byte sequences
#define convert_uint32(x) (((uint32_t)(x)[0]<<24) | ((uint32_t)(x)[1]<<16) | ((uint32_t)(x)[2]<<8) | ((uint32_t)(x)[3]))
#define convert_uint16(x) (((uint16_t)(x)[0]<< 8) | ((uint16_t)(x)[1]))
#define convert_uint8(x)  (((uint8_t) (x)[0]))

// keep the following unsigned, otherwise the sign bit gets messed up
#define convert_int32(x)  (((uint32_t)(x)[0]<<24) | ((uint32_t)(x)[1]<<16) | ((uint32_t)(x)[2]<<8) | ((uint32_t)(x)[3]))
#define convert_int16(x)  (((uint16_t)(x)[0]<< 8) | ((uint16_t)(x)[1]))
#define convert_int8(x)   (((uint8_t) (x)[0]))

#define uint32_byte0(x) (((x) & 0xff000000)>>24)  // MSB
#define uint32_byte1(x) (((x) & 0x00ff0000)>>16)
#define uint32_byte2(x) (((x) & 0x0000ff00)>>8 )
#define uint32_byte3(x) (((x) & 0x000000ff)    )  // LSB

/***************************************************************************************************
 * networking 
 **************************************************************************************************/

// Newer Ethernet shields have a MAC address printed on a sticker on the shield
//byte mac[] = { 
//  0x90, 0xA2, 0xDA, 0x0D, 0x2A, 0x4A }; // this is my ethershield
byte mac[] = { 
  0x04, 0x1E, 0x64, 0x28, 0xD0, 0xE7};  // this is my wiznet w5100 module

char server[] = "192.168.1.100";
int port = 1972;  // 1972 is the default for a FieldTrip buffer

// Initialize the Ethernet client library with the IP address and port of the server 
// that you want to connect to (port 80 is default for HTTP):
EthernetClient client;

/***************************************************************************************************
 * RF12 configuration for 433 Hz 
 **************************************************************************************************/

#define RF12_CS  9
#define RF12_ID  1
#define RF12_GRP 197

/***************************************************************************************************
 * signal processing 
 **************************************************************************************************/

#define sqr(x) ((float)(x)*(float)(x))
#define MAXCHAN 16

float lastval[MAXCHAN];
float filtval[MAXCHAN];
float rmsdecay  = 0.30;     // exponential decay filter for the rms value
float chandecay = 0.10;     // exponential decay filter for the individual channels
float rmsscale = 6000;      // scale the smoothed RMS value between 0-10
float rmsoffset = 800;      // scale the smoothed RMS value between 0-10
float Fvalue = 0;           // filtered value
byte  Bvalue = 0;

/***************************************************************************************************
 * setup
 **************************************************************************************************/

void setup() {

  // setup the signal processing
  for (int i=0; i<MAXCHAN; i++) {
    lastval[i] = 0;
    filtval[i] = 0;
  }

  // Open serial communications and wait for port to open:
  Serial.begin(57600);
  while (!Serial) {
    ; // wait for serial port to connect. Needed for Leonardo and Sparkfun Pro Micro
  }
  Serial.println("[gravity_eeg_network]");

  // start the Ethernet connection:
  while (Ethernet.begin(mac) == 0) {
    Serial.println("Failed to configure Ethernet using DHCP");
    delay(1000);
  }
  Serial.println("Ethernet started");

  // print your local IP address:
  Serial.print("My IP address: ");
  for (byte thisByte = 0; thisByte < 4; thisByte++) {
    // print the value of each byte of the IP address:
    Serial.print(Ethernet.localIP()[thisByte], DEC);
    Serial.print("."); 
  }
  Serial.println();

  // start the RF12 module:
  rf12_set_cs(RF12_CS);
  rf12_initialize(RF12_ID, RF12_433MHZ, RF12_GRP);
  Serial.println("RF12 initialized");

}

/***************************************************************************************************
 * loop
 **************************************************************************************************/

void loop()
{
  while (!client.connected()) {
    // connect to the FieldTrip buffer server
    if (client.connect(server, port)) {
      Serial.println("Connected to FieldTrip buffer");
      break;
    } 
    else {
      Serial.print("Failed to connect to FieldTrip buffer on ");
      Serial.println(server);
      client.stop();
      delay(1000);    
    }
  } // while

  // send the GET_HDR request
  if (client.write(get_hdr, 8)!=8)
    Serial.println("Problem writing request");

  // wait for the response to become available
  unsigned long timeout = millis();
  while (!client.available() && (millis()-timeout)<3000)
    delay(10);

  // get the general messagedef_t section
  for (int i=0; i<8; i++) {
    buf[i] = client.read();
  }

  response.bufsize = convert_uint32(buf+4);
  response.version = convert_uint16(buf+0);
  response.command = convert_uint16(buf+2);

  if ((response.version==VERSION) && (response.command==GET_OK)) {

    // get the headerdef_t section
    for (int i=0; i<response.bufsize; i++) {
      buf[i] = client.read();
    }

    // convert 4 bytes to 32 bit single-precision float
    // note that these do not have to be swapped
    single.bval[3] = buf[12];
    single.bval[2] = buf[13];
    single.bval[1] = buf[14];
    single.bval[0] = buf[15];

    header.nchans    = convert_uint32(buf+0);
    header.nsamples  = convert_uint32(buf+4);
    header.nevents   = convert_uint32(buf+8);
    header.fsample   = single.fval; // see above
    header.data_type = convert_uint32(buf+16);
    header.bufsize   = convert_uint32(buf+20);

    // discard any remaining bytes
    client.flush();

    //    Serial.print("nchans = ");
    //    Serial.print(header.nchans, DEC);
    //    Serial.print("\tnsamples = ");
    //    Serial.print(header.nsamples, DEC);
    //    Serial.print("\tfsample = ");
    //    Serial.print(header.fsample, 2);
    //    Serial.print("\ttype = ");
    //    Serial.print(header.data_type, DEC);
    //    Serial.println();
    //    Serial.println();

    unsigned int nsamples = 1;

    // data selection starts at the desired sample
    get_dat[ 8] = uint32_byte0(header.nsamples-nsamples);
    get_dat[ 9] = uint32_byte1(header.nsamples-nsamples);
    get_dat[10] = uint32_byte2(header.nsamples-nsamples);
    get_dat[11] = uint32_byte3(header.nsamples-nsamples);
    // data selection ends at the last sample (minus 1, since 0 offset)
    get_dat[12] = uint32_byte0(header.nsamples-1);
    get_dat[13] = uint32_byte1(header.nsamples-1);
    get_dat[14] = uint32_byte2(header.nsamples-1);
    get_dat[15] = uint32_byte3(header.nsamples-1);

    // send the GET_DAT request
    if (client.write(get_dat, 16)!=16)
      Serial.println("Problem writing request");

    // wait for the response to become available
    unsigned long timeout = millis();
    while (!client.available() && (millis()-timeout)<2000)
      delay(10);

    // get the general messagedef_t section
    for (int i=0; i<8; i++) {
      buf[i] = client.read();
    }

    response.version = convert_uint16(buf+0);
    response.command = convert_uint16(buf+2);
    response.bufsize = convert_uint32(buf+4);

    //    Serial.print("response.bufsize = ");
    //    Serial.print(response.bufsize, DEC);
    //    Serial.print("\tavailable = ");
    //    Serial.print(client.available(), DEC);
    //    Serial.println();

    // get the datadef_t section
    if ((response.version==VERSION) && (response.command==GET_OK)) {

      for (int i=0; i<response.bufsize; i++) {
        buf[i] = client.read();
      }

      data.nchans    = convert_uint32(buf+0);
      data.nsamples  = convert_uint32(buf+4);
      data.data_type = convert_uint32(buf+8);
      data.bufsize   = convert_uint32(buf+12);

      //      Serial.print("data.bufsize = ");
      //      Serial.print(data.bufsize, DEC);
      //      Serial.println();

      // discard any remaining bytes
      client.flush();

      int offset;
      float value = 0;

      for (unsigned int j=0; j<data.nsamples; j++) {
        for (unsigned int i=0; i<data.nchans; i++) {

          switch (data.data_type) {

          case DATATYPE_CHAR:
          case DATATYPE_UINT8:
            offset = sizeof(datadef_t) + (j*data.nchans+i)*1;
            lastval[i] = convert_uint8(buf+offset);
            break;      

          case DATATYPE_INT8:
            offset = sizeof(datadef_t) + (j*data.nchans+i)*1;
            lastval[i] = convert_int8(buf+offset);
            break;      

          case DATATYPE_UINT16:
            offset = sizeof(datadef_t) + (j*data.nchans+i)*2;
            lastval[i] = convert_uint16(buf+offset);
            break;      

          case DATATYPE_INT16:
            offset = sizeof(datadef_t) + (j*data.nchans+i)*2;
            lastval[i] = convert_int16(buf+offset);
            break;      

          case DATATYPE_UINT32:
            offset = sizeof(datadef_t) + (j*data.nchans+i)*4;
            lastval[i] = convert_uint32(buf+offset);
            break;      

          case DATATYPE_INT32:
            offset = sizeof(datadef_t) + (j*data.nchans+i)*4;
            lastval[i] = convert_int32(buf+offset);
            break;      

          case DATATYPE_FLOAT32:  
            offset = sizeof(datadef_t) + (j*data.nchans+i)*4;
            // convert to 32 bit floating point value            
            // note that these have to be swapped
            single.bval[3] = buf[offset+0];
            single.bval[2] = buf[offset+1];
            single.bval[1] = buf[offset+2];
            single.bval[0] = buf[offset+3];
            lastval[i] = single.fval;
            break;      

          case DATATYPE_FLOAT64:  
          case DATATYPE_UINT64:
          case DATATYPE_INT64:
          default:
            Serial.println("Unsupported data type");
            break;
          } // switch type
        } // loop over channels

        // signal processing
        for (int i=0; i<data.nchans; i++) {
          if (filtval[i]==0)
            filtval[i] = lastval[i];
          else
            filtval[i] = chandecay*(lastval[i]) + (1-chandecay)*filtval[i];
          // compute the sum-of-qsuares
          value += sqr(lastval[i] - filtval[i]);
        }

      } // loop over samples

      // compute the RMS over all channels and samples
      value = sqrt(value/(data.nchans*data.nsamples)); 
      // smooth the RMS with an exponential decay       
      Fvalue = rmsdecay*(value) + (1-rmsdecay)*Fvalue;
      // clip between 0 and 10 and convert to byte
      Bvalue = (byte)(max(min(10*(Fvalue-rmsoffset)/rmsscale, 10), 0));  

      Serial.print(Fvalue);
      Serial.print(" ");
      Serial.print(Bvalue);
      Serial.println();

      unsigned long timeout = millis();

      rf12_sleep(RF12_WAKEUP);
      while(!rf12_canSend() && (millis()-timeout)<1000)
        rf12_recvDone();
      if ((millis()-timeout)<1000) {
        rf12_sendStart(0, &Bvalue, 1);
        rf12_sendWait(0);
      }
      else {
        Serial.println("RF12 timeout, re-initializing");
        rf12_set_cs(RF12_CS);
        rf12_initialize(RF12_ID, RF12_433MHZ, RF12_GRP);
      }
      rf12_sleep(RF12_SLEEP);

    } // got data
    else {
      Serial.println("Error reading data");
      // discard any remaining bytes
      client.flush();
      // dicsonnect from the FieldTrip buffer server
      client.stop();
      delay(1000);
    }

  } // got header
  else {
    Serial.println("Error reading header");
    // discard any remaining bytes
    client.flush();
    // dicsonnect from the FieldTrip buffer server
    client.stop();
    delay(1000);
  }

  // discard any remaining bytes
  client.flush();

  // reading the header and a data packet takes about 32 ms
  // process one packet approximately every 50 ms
  delay(18);
} // loop





