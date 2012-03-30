/*
  Copyright 2012, Donders Center for Cognitive Neuroimaging.
  Based on tia_client_main.cpp, copyright 2010 Graz University of Technology.
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
  [ ] Add (consistent) error handling
  [V] Consider removing thread
  [V] Should we support UDP? -> I think not
  [V] Support "localhost" instead of 127.0.0.1
  [V] Think of support for heterogeneous streams.
  [ ] Add shared-memory access to internal buffer.
  [ ] Ask Christian to read code.
  [ ] Exit cleanly on ctrl-c.
  [ ] Add windows support.
  [ ] Update README with CLI interface.
  [ ] COPYING file for tia2ft, and in libs.
  [ ] Test in TOBI environment.
*/

// STL
#include <iostream>

// Boost
#include <boost/asio.hpp>
#include <boost/thread.hpp>
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

// Forward declarations
int ft_put_data(int ft_buffer_handle, int nchannels, int nsamples, const float
  *chan_samp);
int ft_put_hdr(int ft_buffer_handle, int nchann, int fsample);

bool connect_tia_client(tia::TiAClient &client, const string tia_serv_addr, int
  tia_serv_port);
void start_ft_buffer(int port);
int  sync_meta_info(tia::TiAClient &tia_client, int ft_buffer_handle);
int  forward_packet(tia::DataPacket &packet, int ft_buffer_handle);

int main(int argc, char *argv[])
{
  // Variables for CLI config:
  string tia_host, ft_host;
  int tia_port, ft_port;
  bool serve_ft_buffer, verbose;

  // Since boost is used in TiA's header files, we might use exploit that fact
  // here as well for argument parsing:
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Show help message.")
    ("verbose,v",
     po::value<bool>(&verbose)->default_value(false)->zero_tokens(),
     "Print more info.")
    ("tia-host", po::value<string>(&tia_host)->default_value("localhost"),
     "Set host name of TiA server.")
    ("tia-port", po::value<int>(&tia_port)->default_value(9000),
     "Set port of TiA server.")
    ("fieldtrip-host", po::value<string>(&ft_host)->default_value("localhost"),
     "Set host name of FieldTrip buffer server.")
    ("fieldtrip-port", po::value<int>(&ft_port)->default_value(1972),
     "Set port of FieldTrip buffer server.")
    ("serve-ft-buffer",
     po::value<bool>(&serve_ft_buffer)->default_value(true)->zero_tokens(),
     "Start a new FieldTrip buffer instead of connecting to an existing one.")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  // Show help if requested.
  if (vm.count("help")) {
    cout << desc << endl;
    return 1;
  }

  // Show config:
  cout << "TiA: " << tia_host << ":" << tia_port << endl;
  cout << "FT: " << ft_host << ":" << ft_port << endl;

  // Connect to TiA.
  tia::TiAClient tia_client(true);  // use new-style TiA implementation
  if (!connect_tia_client(tia_client, tia_host, tia_port)) {
    cerr << "Could not connect to TiA server. Closing." << endl;
    cout << desc << endl;
    exit(-1);
  }
  cout << "Connected to TiA server." << endl;
  tia_client.startReceiving(false);  // use TCP


  // Connect to FT buffer, and start one if requested:
  int ft_buffer_handle = 0;  // zero indicates shared-memory communication.
  if (serve_ft_buffer) {
    boost::thread ft_buffer_thread(&start_ft_buffer, ft_port);
    cout << "Started FieldTrip buffer on port " << ft_port << "." << endl;
    // TODO: implement tearing down the FT-buffer when supported in libbuffer.
  } else {
    ft_buffer_handle = open_connection(ft_host.c_str(), ft_port);
    if(ft_buffer_handle <= 0) {
      cerr << "Could not connect to FieldTrip buffer. Closing." << endl;
      cout << desc << endl;
      exit(-1);
    }
  }

  // Sync information on data stream:
  sync_meta_info(tia_client, ft_buffer_handle);

  // Main loop
  tia::DataPacket *packet = tia_client.getEmptyDataPacket();
  while (true) {
    cout << "." << flush;
    tia_client.getDataPacket(*packet);

    if (verbose) {
      cout << "#: " << packet->getConnectionPacketNr() << endl;
      cout << "time: " << packet->getTimestamp() << endl;
      cout << "nchan: " << packet->getNrOfChannels()[0] << endl;
      cout << "samp/chan: " << packet->getNrSamplesPerChannel()[0] << endl;
    }

    forward_packet(*packet, ft_buffer_handle);
  }

  try {
    tia_client.stopReceiving();
  } catch (std::exception &e) {
    cerr << "Stop Receiving failed -- Error:" << "--> " << e.what() << endl;
  }

  return(0);
}


/* Connect to TiA server. */
bool connect_tia_client(tia::TiAClient &client, const string tia_serv_addr,
                        int tia_serv_port)
{
  try {
    client.connect(tia_serv_addr == "localhost" ? "127.0.0.1" : tia_serv_addr, 
      tia_serv_port);
    client.requestConfig();
  } catch(std::exception &e) {
    cerr << e.what() << endl;
    return 1;
  }
}


/* Starts FieldTrip buffer server on this host. */
void start_ft_buffer(int port)
{
  host_t host;
  host.port = port;
  check_datatypes();  // sanity check for sizes of different data types.
  tcpserver((void *)(&host));  // start server
}


/* Read meta information (channel configuration etc.) from TiA, and write this
 * to a FieldTrip buffer. So far, it only *prints* the meta info. */
int sync_meta_info(tia::TiAClient &tia_client, int ft_buffer_handle)
{
  // Request TiA config:
  try {
    tia_client.requestConfig();
  } catch (std::exception &e) {
    cerr << "Requesting TiA config failed -- Error:" << "--> "
         << e.what() << endl;
    return 1;
  }

  tia::SignalInfo sigInfo = tia_client.config().signal_info;

  // Print some signal statistics
  float fsample = sigInfo.masterSamplingRate();
  cout << "Detected the following meta information:" << endl;
  cout << "Sampling rate: " << fsample << endl;
  cout << "Block size: " << sigInfo.masterBlockSize() << endl;

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

  cout << "Sending to FieldTrip buffer..." << endl;
  return ft_put_hdr(ft_buffer_handle, nchann, fsample);
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
 *
 * Note that even if different signal streams have the same sampling rate, TiA
 * does not guarantee that an equal amount of samples per signal stream is
 * present in a data packet.
 */
int forward_packet(tia::DataPacket &packet, int ft_buffer_handle)
{
  if (packet.getNrOfSignalTypes() != 1) {
    cerr << "Heterogeneous signal streams are not yet supported :/." << endl;
    exit(2);
  }

  /* Vector tia_raw contains the data. The values are multiplexed as follows:
   * [channel 1 sample 1...n, channel 2 sample 1..n, ...]
   */
  std::vector<double> tia_raw = packet.getData();
  int nchannels = packet.getNrOfChannels()[0];
  int nsamples = packet.getNrSamplesPerChannel()[0];

  // Convert to FT-buffer byte-order.
  std::vector<float> ft_raw(tia_raw.size(), 0);
  for (int ci=0; ci < nchannels; ++ci)
    for (int si=0; si < nsamples; ++si)
      ft_raw[ci + si * nchannels] = tia_raw[si + ci * nsamples];

  return ft_put_data(ft_buffer_handle, nchannels, nsamples, &ft_raw[0]);
}


/* Convenience function to call put data (PUT_DAT) in a FT buffer */
int ft_put_hdr(int ft_buffer_handle, int nchann, int fsample)
{
  // Build data-structures for sending header information.
  headerdef_t hdr = {0};
  hdr.nchans = nchann;
  hdr.fsample = fsample;
  hdr.data_type = DATATYPE_FLOAT32;

  messagedef_t req_hdr = {0};
  req_hdr.version = VERSION;  // VERSION could use a prefix.
  req_hdr.command = PUT_HDR;
  req_hdr.bufsize = sizeof(req_hdr) + sizeof(headerdef_t);

  message_t req = {0};
  req.def = &req_hdr;
  req.buf = malloc(req_hdr.bufsize);
  memset(req.buf, 0, req_hdr.bufsize);
  memcpy(req.buf, &hdr, sizeof(hdr));
  // TODO add FT_CHUNK_CHANNEL_NAMES

  message_t *response = NULL;  // client request expects a *double* pointer
  int status = clientrequest(ft_buffer_handle, &req, &response);
  cout << "clientreq returned: " << status << endl;
  assert(response->def->command == PUT_OK);

  // Note: if PUT_HDR fails, it still can modify the HDR information.

  free(req.buf);
  free(response->buf);
  free(response->def);
  free(response);
}


/* Convenience wrapper to perform a PUT_DAT on ft_buffer */
int ft_put_data(int ft_buffer_handle, int nchannels, int nsamples,
                const float *chan_samp)
{
  int status = 0;
  // Create descriptor of raw data for FT-buffer:
  datadef_t data_hdr = {0};
  data_hdr.nchans = nchannels;
  data_hdr.nsamples = nsamples;
  data_hdr.data_type = DATATYPE_FLOAT32;
  data_hdr.bufsize = nchannels * nsamples * sizeof(float);

  // Create packet header for FT-buffer request
  messagedef_t req_hdr = {0};
  req_hdr.version = VERSION;
  req_hdr.command = PUT_DAT;
  req_hdr.bufsize = sizeof(data_hdr) + data_hdr.bufsize;

  // Compose final FT-buffer request
  message_t req = {0};
  req.def = &req_hdr;
  req.buf = malloc(req_hdr.bufsize);
  memset(req.buf, 0, req_hdr.bufsize);
  memcpy(req.buf, &data_hdr, sizeof(data_hdr));
  memcpy((char *) req.buf + sizeof(data_hdr), chan_samp, data_hdr.bufsize);

  message_t *response = NULL;
  if (!clientrequest(ft_buffer_handle, &req, &response))
    status = -1;
  if (response->def->command != PUT_OK)
    status = -2;

  free(req.buf);
  free(response->buf);
  free(response->def);
  free(response);

  return status;
}

