/*  Copyright 2010 Graz University of Technology
    Contact: SignalServer@tobi-project.org
*/

// STL
#include <iostream>

// Boost
#include <boost/asio.hpp>

// TiA
#include "tia/data_packet_interface.h"
#include "tia/tia_client.h"
#include "tia/ssconfig.h"
#include "tia/defines.h"

// FieldTrip buffer
extern "C" { // TODO: is this needed?
#include "buffer.h"
}

using namespace std;

bool connect_tia_client(tia::TiAClient &client, const string srv_addr, 
  int srv_port) {
  // connect to TiA server
  try {
    client.connect(srv_addr, srv_port);
    client.requestConfig();
  } 
  catch(std::exception &e) {
    cerr << e.what() << endl;
    return 1;
  }
}


/* Starts FieldTrip buffer server on this host. */
void start_ft_buffer(int port) {
  // TODO: start in separate thread
  host_t host;
  host.port = port;
  check_datatypes();  // sanity check for sizes of different data types.
  tcpserver((void *)(&host));  // start server
}


/* Read meta information (channel configuration etc.) from TiA, and write this
 * to a FieldTrip buffer. So far, it only *prints* the meta info. */
int sync_meta_info(tia::TiAClient &tia_client, void *ft_buffer) {
  // Request TiA config:
  try {
    tia_client.requestConfig();
  } 
  catch (std::exception &e) {
    cerr << "Requesting config failed -- Error:" << "--> " << e.what() << endl;
    return 1;
  }

  tia::SignalInfo sigInfo = tia_client.config().signal_info;

  // Print some signal statistics
  cout << "Detected the following meta information:" << endl;
  cout << "Sampling rate: " << sigInfo.masterSamplingRate() << endl;
  cout << "Block size: " << sigInfo.masterBlockSize() << endl;

  tia::SignalInfo::SignalMap::iterator i;
  for(i = sigInfo.signals().begin(); i != sigInfo.signals().end(); ++i) {
    // i contains (string, Signal) pairs
    cout << i->first << endl;

    tia::Signal signal (i-> second);
    cout << "Signal (modality) fs: " << signal.samplingRate() << endl;

    std::vector<tia::Channel> channels = signal.channels();
    cout << "Detected " << channels.size() << " channels:" << endl;
    for(int i = 0; i < channels.size(); ++i) {
      cout << channels[i].id() << (i < channels.size() - 1 ? ", " : "");
    }
  }
  cout << "." << endl;  // end enumeration.

  // TODO: Send info to FT-buffer
}

int main(int argc, const char *argv[])
{
  /*
    Outline:
    1) Parse command line output,
    2) connect to TiA server,
    3) connect to FT buffer,
    4) read header, put header
    5) start streaming, put data in buffer

    TODO:
    [ ] add (consistent) error handling
    [V] consider removing thread
    [ ] should we support UDP? -> I think not
  */

  string srv_addr = "127.0.0.1";
  boost::uint16_t srv_port = 9000;
  bool host_ft_buffer = false;

  // parse command line
  if(argc == 1) {
    cout << "Using default server " << srv_addr << ":" << srv_port << endl;
  } else if(argc == 2) {
    string param(argv[1]);
  } else if(argc == 3) {
    srv_addr = argv[1];
    stringstream conv(argv[2]);
    conv >> srv_port;
    cout << "Using server " << srv_addr << ":" << srv_port << endl;
  } else if(argc == 4) {
    string param(argv[1]);
    srv_addr = argv[2];
    stringstream conv(argv[3]);
    conv >> srv_port;
    cout << "Using server " << srv_addr << ":" << srv_port << endl;
  } else {
    cout << "Wrong number of arguments given: " << argc-1 << endl;
    cout << " - Usage: " << argv[0] << "  signalserver-ip   port" << endl;
    return(-1);
  }

  // Initialize connections
  tia::TiAClient client(true);  // use new-style TiA implementation
  connect_tia_client(client, srv_addr, srv_port);
  cout << "Connected to TiA server." << endl;

  if (host_ft_buffer)
    start_ft_buffer(1972);

  // Start receiving:
  client.startReceiving(false);  // use TCP
  tia::DataPacket *packet = client.getEmptyDataPacket();

  // TODO connect to FT buffer

  sync_meta_info(client, NULL);

  // Main loop
  while (true) {
    cout << "Getting packet..." << flush;
    client.getDataPacket(*packet);

    cout << "#: " << packet->getConnectionPacketNr() << endl;
    cout << "Flags: " << packet->getFlags() << endl;
    cout << "time: " << packet->getTimestamp() << endl;

    assert(packet->getNrOfSignalTypes() == 1);  // no support for mixed yet
    cout << "nchan: " << packet->getNrOfChannels()[0] << endl;
    cout << "nsamp: " << packet->getNrOfSamples()[0] << endl;
    cout << "samp/chan: " << packet->getNrSamplesPerChannel()[0] << endl;

    std::vector<double> const payload(packet->getData());
    cout << "|data|: " << payload.size() << endl;
  }

  try {
    client.stopReceiving();
  } catch (std::exception &e) {
    cerr << "Stop Receiving failed -- Error:" << "--> " << e.what() << endl;
  }

  return(0);
}
