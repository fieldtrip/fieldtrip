/*  Copyright 2010 Graz University of Technology
    Contact: SignalServer@tobi-project.org
*/

// STL
#include <iostream>
#include <algorithm>

// Boost
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/thread.hpp>
#include <boost/cstdint.hpp>

// local
#include "tia/data_packet_interface.h"
#include "tia/tia_client.h"
#include "tia/ssconfig.h"
#include "tia/defines.h"

using namespace std;

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
    [ ] add error handling
    [ ] consider removing thread
    [ ] should we support UDP? -> I think not
  */

  string srv_addr = "127.0.0.1";
  boost::uint16_t srv_port = 9000;

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

  tia::TiAClient client(true);

  // connect to TiA server
  try {
    client.connect(srv_addr, srv_port);
    client.requestConfig();
    cout << "Connected to TiA server." << endl;
  } catch(std::exception &e) {
    cerr << e.what() << endl;
    return 1;
  }

  // TODO: connect to FieldTrip buffer

#ifdef WIN32
  SetPriorityClass(GetCurrentProcess(), ABOVE_NORMAL_PRIORITY_CLASS);
  SetPriorityClass(reader_thread.native_handle(), REALTIME_PRIORITY_CLASS);
  SetThreadPriority(reader_thread.native_handle(), THREAD_PRIORITY_TIME_CRITICAL );
#endif

  // Request TiA config:
  try {
    client.requestConfig();
    tia::SignalInfo sigInfo = client.config().signal_info;

    // Print some signal statistics
    cout << "Sampling rate: " << sigInfo.masterSamplingRate() << endl;
    cout << "Block size: " << sigInfo.masterBlockSize() << endl;

    tia::SignalInfo::SignalMap sigMap = sigInfo.signals();
    tia::SignalInfo::SignalMap::iterator i;
    for(i = sigMap.begin(); i != sigMap.end(); ++i) {
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
    cout << "." << endl;
  } catch (std::exception &e) {
    cerr << "Requesting config failed -- Error:" << "--> " << e.what() << endl;
  }

  // Start receiving:
  client.startReceiving(false);  // use TCP
  tia::DataPacket *packet = client.getEmptyDataPacket();

  // Main loop
  while (true) {
    cout << "Getting packet..." << flush;
    client.getDataPacket(*packet);
    cout << "yay!" << endl;


    // TODO learn to use and inspect packet.
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
