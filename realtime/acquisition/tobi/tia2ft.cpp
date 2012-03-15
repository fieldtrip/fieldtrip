/* 
  Copyright 2012, Donders institute of Cognitive Neuroimaging.
  Based on tia_clien_main.cpp, copyright 2010 Graz University of Technology.
  License: GPL or BSD.
*/

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
  [ ] support "localhost" instead of 127.0.0.1
  [V] think of support for heterogeneous streams.
  [ ] think of support for heterogeneous streams.
  [ ] ask Christian to read code.[
  [ ] test in TOBI environment.
*/

// STL
#include <iostream>

// Boost
#include <boost/asio.hpp>
#include <boost/program_options.hpp>

// TiA
#include "tia/data_packet_interface.h"
#include "tia/tia_client.h"
#include "tia/ssconfig.h"
#include "tia/defines.h"

// FieldTrip buffer
extern "C" {
#include "buffer.h"
}

using namespace std;
namespace po = boost::program_options;

bool connect_tia_client(tia::TiAClient &client, const string tia_serv_addr, 
  int tia_serv_port) {
  // connect to TiA server
  try {
    client.connect(tia_serv_addr, tia_serv_port);
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
int sync_meta_info(tia::TiAClient &tia_client, int ft_buffer) {
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
  float fsample = sigInfo.masterSamplingRate();
  cout << "Detected the following meta information:" << endl;
  cout << "Sampling rate: " << fsample << endl;
  cout << "Block size: " << sigInfo.masterBlockSize() << endl;

  /*
  int nchann = 0;
  tia::SignalInfo::SignalMap::iterator i;
  for(i = sigInfo.signals().begin(); i != sigInfo.signals().end(); ++i) {
    // i contains (string, Signal) pairs
    cout << i->first << endl;

    tia::Signal signal (i->second);
    cout << "Signal (modality) fs: " << signal.samplingRate() << endl;

    std::vector<tia::Channel> channels = signal.channels();
    cout << "Detected " << channels.size() << " channels:" << endl;
    nchann += channels.size();
    for(int i = 0; i < channels.size(); ++i) {
      cout << channels[i].id() << (i < channels.size() - 1 ? ", " : "");
    }
  }
  cout << "." << endl;  // end enumeration.
  */

  tia::Signal signal = sigInfo.signals().begin()->second;
  int nchann = signal.channels().size();

  /* Attempt at shorter extraction of channel labels
  stringstream sens_lab;
  for(int i = 0; i < nchann; ++i) {
    sens_lab 
      //<< signal.channels[i].id() 
      << "\n";
  }
  
  cout << sens_lab.c_str() << endl;
  */

  // TODO: Send info to FT-buffer
  cout << "Sending to FieldTrip buffer..." << endl;
  
  // Build data-structures for sending header information.

  headerdef_t hdr = {0};
  hdr.nchans = nchann;
  hdr.fsample = fsample;
  hdr.data_type = DATATYPE_FLOAT32;

  messagedef_t req_hdr = {0};
  req_hdr.version = VERSION;  // VERSION could use a prefix.
  req_hdr.command = PUT_HDR;

  message_t req = {0};
  req.def = &req_hdr;
  req.buf = &hdr;
  req_hdr.bufsize = sizeof(headerdef_t);  // update with buffer size
  // add FT_CHUNK_CHANNEL_NAMES

  message_t *response = NULL;  // client request expects a *double* pointer
  int status = clientrequest(ft_buffer, &req, &response);
  cout << "clientreq returned: " << status << endl;
  assert(response->def->bufsize == 0);  // unexpected payload in response
  // TODO: free response (check for mem-leaks)
}


/* Take a TiA packet, and append it in the FieldTrip buffer. Apparently, the
 * packet *cannot be const*. 
 * 
 * Note that TiA support streaming of heterogeneous data streams, that can have
 * different sampling rates and block-sizes. Since the FT-buffer only supports
 * homogeneous data streams and there is no demand yet, only homogeneous
 * data streams are supported.
 * 
 * Support for heterogeneous data streams can be implemented by using a
 * sampling rate that an integer multiple of the modality's sampling rates, and
 * re-sampling the signals accordingly.
 */
void forward_packet(tia::DataPacket &packet, int ft_buffer)
{
  if (packet.getNrOfSignalTypes() != 1){
    cerr << "Heterogeneous signal streams are not yet supported :/." << endl;
    exit(2);
  }

  std::vector<float> const payload(
    packet.getData().begin(), packet.getData().end());  // cast to float

  cout << "|data|: " << payload.size() << endl;

  /* Vector payload now contains the data. The values are multiplexed as
   * follows: 
   * [Ch1_samp1, Ch2_samp1, ... Chn_samp1, Ch1_samp2, Chan2_samp2, ...].
   */

  // Add raw data to FT-buffer
  datadef_t data_hdr;
  memset(&data_hdr, 0, sizeof(datadef_t));
  data_hdr.nchans = packet.getNrOfChannels()[0];
  data_hdr.nsamples = packet.getNrOfSamples()[0];
  data_hdr.data_type = DATATYPE_FLOAT32; // FIXME: already stated in PUT_HDR?
  data_hdr.bufsize = payload.size() * sizeof(float);  // Redundant?
  cout << "Creating packet of " << data_hdr.nchans << "x" 
    << data_hdr.nsamples << "@" << sizeof(float) << " bytes." << endl;

  data_t data = {0};
  data.def = &data_hdr;
  data.buf = (void *) &payload[0]; // TODO: cast & find correct accessor.
  cout << "ft packet payload: " << payload[0] << endl;

  messagedef_t req_hdr = {0};
  req_hdr.version = VERSION;
  req_hdr.command = PUT_DAT;

  message_t req = {0};
  req.def = &req_hdr;
  req.buf = &data;
  req_hdr.bufsize = sizeof(data) + data.def->bufsize; 
  cout << "req_hdr.bufsize = " << req_hdr.bufsize << endl;

  message_t *response = NULL; 
  int status = clientrequest(ft_buffer, &req, &response);
  cout << "clientreq returned: " << status << endl;
  assert(status == 0);
  assert(response->def->bufsize == 0);  // unexpected payload in response
  // TODO: free response (check for mem-leaks)
}

int main(int argc, char *argv[])
{

  // Since boost is used in TiA's header files, we might use exploit that fact
  // here as well for argument parsing:

  string tia_host, ft_host;
  int tia_port, ft_port;
  bool serve_ft_buffer;

  // parse command line
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Show help message.")
    ("tia-host", po::value<string>(&tia_host)->default_value("127.0.0.1"), 
      "Set host name of TiA server.")
    ("tia-port", po::value<int>(&tia_port)->default_value(9000), 
      "Set port of TiA server.")
    ("fieldtrip-host", po::value<string>(&ft_host)->default_value("localhost"), 
      "Set host name of FieldTrip buffer server.")
    ("fieldtrip-port", po::value<int>(&ft_port)->default_value(1972), 
      "Set port of FieldTrip buffer server.")
    ("serve-ft-buffer", po::value<bool>(&serve_ft_buffer)->default_value(false),
      "Start a new FieldTrip buffer instead of connecting to an existing one.")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << desc << endl;
    return 1;
  }

  cout << "tia: " << tia_host << ":" << tia_port << endl;
  cout << "ft: " << ft_host << ":" << ft_port << endl;

  // initialize connections
  tia::TiAClient tia_client(true);  // use new-style TiA implementation
  connect_tia_client(tia_client, tia_host, tia_port);
  cout << "Connected to TiA server." << endl;

  if (serve_ft_buffer)
    start_ft_buffer(1972);

  // Start receiving:
  tia_client.startReceiving(false);  // use TCP

  // TODO connect to FT buffer
  int ft_buffer = open_connection(ft_host.c_str(), ft_port);
  if(ft_buffer <= 0)
  {
    cerr << "Could not connect to FieldTrip buffer. Closing." << endl;
    exit(-1);
  }

  sync_meta_info(tia_client, ft_buffer);

  // Main loop
  tia::DataPacket *packet = tia_client.getEmptyDataPacket();
  while (true) {
    cout << "Getting packet..." << flush;
    tia_client.getDataPacket(*packet);

    cout << "#: " << packet->getConnectionPacketNr() << endl;
    cout << "Flags: " << packet->getFlags() << endl;
    cout << "time: " << packet->getTimestamp() << endl;
    cout << "nchan: " << packet->getNrOfChannels()[0] << endl;
    cout << "nsamp: " << packet->getNrOfSamples()[0] << endl;
    cout << "samp/chan: " << packet->getNrSamplesPerChannel()[0] << endl;

    forward_packet(*packet, ft_buffer);
  }

  try {
    tia_client.stopReceiving();
  } catch (std::exception &e) {
    cerr << "Stop Receiving failed -- Error:" << "--> " << e.what() << endl;
  }

  return(0);
}
