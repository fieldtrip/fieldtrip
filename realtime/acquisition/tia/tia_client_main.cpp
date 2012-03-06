/*
    This file is part of the TOBI SignalServer.

    Commercial Usage
    Licensees holding valid Graz University of Technology Commercial
    licenses may use this file in accordance with the Graz University
    of Technology Commercial License Agreement provided with the
    Software or, alternatively, in accordance with the terms contained in
    a written agreement between you and Graz University of Technology.

    --------------------------------------------------

    GNU General Public License Usage
    Alternatively, this file may be used under the terms of the GNU
    General Public License version 3.0 as published by the Free Software
    Foundation and appearing in the file gpl.txt included in the
    packaging of this file.  Please review the following information to
    ensure the GNU General Public License version 3.0 requirements will be
    met: http://www.gnu.org/copyleft/gpl.html.

    In case of GNU General Public License Usage ,the TOBI SignalServer
    is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the TOBI SignalServer. If not, see <http://www.gnu.org/licenses/>.

    Copyright 2010 Graz University of Technology
    Contact: SignalServer@tobi-project.org
*/

/**
* @file ssclient_main.cpp
* @brief This file includes a demo implementation of a very simple console TiA client.
**/

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

//-----------------------------------------------------------------------------

class TiAClientDataReader
{
  public:
    TiAClientDataReader(tia::TiAClient& client, boost::mutex& mutex, boost::condition_variable& cond) :
        client_(client), packet_(0), mutex_(mutex), cond_(cond),
        running_(1),
        timestamp_(boost::posix_time::microsec_clock::local_time()),
        t_min_total_ (10, 0, 0),
        t_max_total_ (0, 0, 0),
        t_min_last_ (10, 0, 0),
        t_max_last_ (0, 0, 0),
        t_var_(0)
    {
      packet_ = client_.getEmptyDataPacket();
    }

    void stop()
    {
      running_ = 0;
    }

    void readData()
    {
      unsigned int counter = 0;

      boost::uint64_t sample_nr = 0;
      boost::uint64_t packet_nr = 0;

      boost::uint64_t sample_nr_old = 0;
      boost::uint64_t packet_nr_old = 0;

      while(running_)
      {
        {

          boost::unique_lock<boost::mutex> lock(mutex_);

          vector<double> v;

          if (!client_.receiving())
          {
            cond_.wait(lock);
          }

          if(running_ && client_.receiving())
          {
            try {

              std::cout << " ---- " << std::endl;

              client_.getDataPacket(*packet_);

              std::cout << "NR SigTypes: " << packet_->getNrOfSignalTypes() << std::endl;
              std::cout << "Flags: " << packet_->getFlags() << std::endl;
              std::cout << "ID: " << packet_->getPacketID() << std::endl;

              std::cout << " *** " << std::endl << std::endl;

              if(packet_->hasFlag(SIG_MOUSE))
              {
                v = packet_->getSingleDataBlock(SIG_MOUSE);
                cout << "mousepos: "<< packet_->getSingleDataBlock(SIG_MOUSE)[1]<<","<<packet_->getSingleDataBlock(SIG_MOUSE)[2];
              }

              if(packet_->hasFlag(SIG_MBUTTON))
              {
                v = packet_->getSingleDataBlock(SIG_MBUTTON);
                cout<<" ... buttons: "<<packet_->getSingleDataBlock(SIG_MBUTTON)[1]<<packet_->getSingleDataBlock(SIG_MBUTTON)[2]<<packet_->getSingleDataBlock(SIG_MBUTTON)[3]<<endl;
              }
              if( ((client_.config().signal_info.masterSamplingRate()/client_.config().signal_info.masterBlockSize()) < 1) ||
                (counter%(
                (client_.config().signal_info.masterSamplingRate()/client_.config().signal_info.masterBlockSize()) *2 ) == 0) )
              {
                std::cout << " -- Nr:" << packet_->getPacketID() << std::endl;
                //                client_.requestConfig();
                //                cerr << " -- got config ... "  << endl;
              }

              sample_nr_old = sample_nr;
              packet_nr_old = packet_nr;

              sample_nr = packet_->getPacketID();
              packet_nr = packet_->getConnectionPacketNr();

              if( (sample_nr - sample_nr_old) > 1)
              {
                cerr << "SampleNr difference: " << (sample_nr - sample_nr_old);
                cerr << " -- @Sample: " << sample_nr << endl;
              }

              if( (packet_nr - packet_nr_old) > 1)
              {
                cerr << "PacketNr difference: " << (packet_nr - packet_nr_old);
                cerr << " -- @Packet: " << packet_nr << endl;
              }

            }
            catch (std::exception& e)
            {

              cerr << "*** " << e.what() << endl;
              break;
            }

            #ifdef TIMING_TEST
              timestamp_ = boost::posix_time::microsec_clock::local_time();
              diff_ = timestamp_ - packet_->getTimestamp();

              t_diffs_.push_back (diff_);
              t_min_total_ = min (t_min_total_, diff_);
              t_max_total_ = max (t_max_total_, diff_);
              t_min_last_ = min (t_min_last_, diff_);
              t_max_last_ = max (t_max_last_, diff_);

              t_mean_ = (t_mean_ + diff_)/2;
              t_var_  = (t_var_ +
              ( (diff_.total_microseconds() - t_mean_.total_microseconds() )*
                (diff_.total_microseconds() - t_mean_.total_microseconds() )  ) )/2;
            #endif
          }
          else
            break;

          counter++;
          #ifdef TIMING_TEST
            if( ((client_.config().signal_info.masterSamplingRate()/client_.config().signal_info.masterBlockSize()) < 1) ||
              (counter%(
              (client_.config().signal_info.masterSamplingRate()/client_.config().signal_info.masterBlockSize()) *2 ) == 0) )
            {
              sort (t_diffs_.begin(), t_diffs_.end());
              cout << "Packet Nr.: " << counter << ";  ";
              cout << "Timing (microsecs) -- mean: " << t_mean_.total_microseconds() << ", ";
              cout << "variance: " << t_var_;
              cout << ", min: " << t_min_last_.total_microseconds() << " (total: "<<  t_min_total_.total_microseconds() <<"), ";
              cout << "max: "<< t_max_last_.total_microseconds() << " (total: "<< t_max_total_.total_microseconds() << "), ";
              cout << "median: " << t_diffs_[t_diffs_.size() / 2].total_microseconds () << endl;
              t_diffs_.clear();
              t_min_last_ = boost::posix_time::time_duration (10, 0, 0);
              t_max_last_ = boost::posix_time::time_duration (0, 0, 0);
            }
          #endif
        }
      }
    }
  private:
    tia::TiAClient&                   client_;
    tia::DataPacket*                  packet_;
    boost::mutex&               mutex_;
    boost::condition_variable&  cond_;
    bool                        running_;
    boost::posix_time::ptime timestamp_;
    boost::posix_time::time_duration diff_;
    boost::posix_time::time_duration t_mean_;
    boost::posix_time::time_duration t_min_total_;
    boost::posix_time::time_duration t_max_total_;
    boost::posix_time::time_duration t_min_last_;
    boost::posix_time::time_duration t_max_last_;
    vector<boost::posix_time::time_duration> t_diffs_;
    boost::int64_t t_var_;
//     boost::posix_time::time_duration t_var_;
};


//-----------------------------------------------------------------------------

int main(int argc, const char* argv[])
{
  string   srv_addr = "127.0.0.1";
  boost::uint16_t srv_port = 9000;

  bool use_new_tia = true;

  if(argc == 1)
  {
    cout << "Using default server " << srv_addr << ":" << srv_port << endl;
  }
  else if(argc == 2)
  {
    string param (argv[1]);
    if (param == "-o")
      use_new_tia = false;
  }
  else if(argc == 3)
  {
    srv_addr = argv[1];
    stringstream conv(argv[2]);
    conv >> srv_port;
    cout << "Using server " << srv_addr << ":" << srv_port << endl;
  }
  else if(argc == 4)
  {
    string param (argv[1]);
    if (param == "-o")
      use_new_tia = false;

    srv_addr = argv[2];
    stringstream conv(argv[3]);
    conv >> srv_port;
    cout << "Using server " << srv_addr << ":" << srv_port << endl;
  }
  else
  {
    cout << "Wrong number of arguments given: " << argc-1 << endl;
    cout << " - Usage: " << argv[0] << "  signalserver-ip   port" << endl;
    return(-1);
  }

  tia::TiAClient client(use_new_tia);

  try
  {
    client.connect(srv_addr, srv_port);
    client.requestConfig();
  }
  catch(std::exception& e)
  {
    cerr << e.what() << endl;
    return 1;
  }

  boost::mutex mutex;
  boost::condition_variable cond;

  TiAClientDataReader reader(client, mutex, cond);
  boost::thread reader_thread(boost::bind(&TiAClientDataReader::readData, &reader));

    #ifdef WIN32
      SetPriorityClass(GetCurrentProcess(), ABOVE_NORMAL_PRIORITY_CLASS);
      SetPriorityClass(reader_thread.native_handle(), REALTIME_PRIORITY_CLASS);
      SetThreadPriority(reader_thread.native_handle(), THREAD_PRIORITY_TIME_CRITICAL );
    #endif

  map<string, string> known_commands;
  known_commands.insert(make_pair("config",   "Requests the config from server."));
  known_commands.insert(make_pair("starttcp", "Starts the data transmission using tcp"));
  known_commands.insert(make_pair("startudp", "Starts the data transmission using udp broadcast"));
  known_commands.insert(make_pair("stop",     "Stops the data transmission"));
  known_commands.insert(make_pair("q",        "Quits the client."));
  known_commands.insert(make_pair("help",     "Prints this help text."));
  known_commands.insert(make_pair("state", "Get the data receiving state of the client"));

  string command;
  cout << endl << ">>";

  while(cin >> command)
  {
    {
      map<string, string>::const_iterator it = known_commands.find(command);
      if (it == known_commands.end())
      {
        cout << "Type 'help' to get description of supported commands" << endl;
        cout << endl << ">>";
        continue;
      }
    }

//    boost::unique_lock<boost::mutex> lock(mutex);

    if(command == "q")
    {
      break;
    }

    if (command == "config")
    {
      try {
        client.requestConfig();
      }
      catch (std::exception& e)
      {
        cerr << "Requesting config failed -- Error:"
             << "--> " << e.what() << endl;
      }
    }
    else if (command == "state")
    {
        cout << "Client state: ";
        if (client.receiving())
            cout << "receiving data" << endl;
        else
            cout << "not receiving data" << endl;
    }
    else if (command == "starttcp" || command == "startudp")
    {
      bool udp =  command == "startudp";

      cout << "Start Receiving ..." << endl;

      try {
        client.startReceiving(udp);
        if (!client.receiving())
            cerr << "Client did not block until data receiving..." << endl;
      }
      catch (std::exception& e)
      {
        cerr << "Start Receiving failed -- Error:" << endl
             << "--> " << e.what() << endl;
      }

      cond.notify_one();
    }
    else if (command == "stop")
    {
      cout << "Stop Receiving ..." << endl;

      //boost::unique_lock<boost::mutex> lock(mutex);

      try {
        client.stopReceiving();
      }
      catch (std::exception& e)
      {
        cerr << "Stop Receiving failed -- Error:" << endl
             << "--> " << e.what() << endl;
      }
    }
    else if (command == "help")
    {
      map<string, string>::const_iterator it = known_commands.begin();
      map<string, string>::const_iterator end = known_commands.end();
      for (; it != end; ++it)
      {
        const string& command_name = (*it).first;
        string spacing;
        spacing.assign(20 - command_name.size(), '.');
        cout << (*it).first << spacing << (*it).second << endl;
      }
    }

    cond.notify_one();
    cout << endl << ">>";
  }

  cout << "Terminating ..." << endl;
  reader.stop();
  cond.notify_all();
  reader_thread.join();

  try {
    client.stopReceiving();
  }
  catch (std::exception& e)
  {
    cerr << "Stop Receiving failed -- Error:"
         << "--> " << e.what() << endl;
  }

 return(0);
}

//-----------------------------------------------------------------------------
