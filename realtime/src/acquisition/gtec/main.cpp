//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#include <iostream>
#include <string>
#include <iomanip>

#include "EstablishNetworkConnection.hpp"
#include "SetupXMLMessage.hpp"
#include "Transceiver.hpp"
#include "ReadMeasurementData.hpp"
#include "EscapeXML.hpp"
#include "XMLSTLConversions.hpp"
#include "SetupgUSBampXMLConfig.hpp"
#include "SetupgHIampXMLConfig.hpp"
#include "SetupgNautilusXMLConfig.hpp"
#include "CheckServerReply.hpp"
#include "Device.hpp"

#include "interface.h"

#define EXPECTED_ARGC 7

//forward declarations (implementation follows main)
std::string ToString(int v);
bool CheckDeviceIndex(int index, const size_t& max_index);
bool CheckSampleRate(int sample_rate, const Device& device);
bool CheckPositive(int number, const size_t& dummy);
template<typename T1, typename T2> void PromptUser(std::ostream* outstream, std::istream* instream, std::string message, std::string error_message, T1* result, bool (*constrain_function)(T1, const T2&), T2 constrain = 0);
void PrintDevices(std::ostream* outstream, const std::vector<Device>& device_list);

// FieldTrip buffer
int ft_server;

//main entry point
int main(int argc, char* argv[])
{
	int ret_main = 0;
	std::cout << "g.NEEDaccess Interface to FieldTrip buffer" << std::endl << std::endl;

	// If there are too less command line arguments, display the usage of the
	// program and return
	//------------------------------------------------------------------------------
	if (argc < EXPECTED_ARGC)
	{
		std::cerr << "Error: Expected " << EXPECTED_ARGC - 1 << " arguments, but got " << argc - 1 << std::endl;
		std::cerr << "Usage: ProgramName HostIp HostPort LocalIp LocalPort BufferIp BufferPort" << std::endl;
		std::cerr << std::endl << "Please press <return> to exit" << std::endl;
		std::cin.get();
		return -1;
	}

	// Parse command line arguments from console to retrieve the program's
	// configuration
	//------------------------------------------------------------------------------
	std::string host_ip(argv[1]);
	int host_port = atoi(argv[2]);
	std::string local_ip(argv[3]);
	int local_port = atoi(argv[4]);
	std::string buffer_ip(argv[5]);
	int buffer_port = atoi(argv[6]);

	// Display the given parameters
	//------------------------------------------------------------------------------
	std::cout << "HostIp     = " << host_ip << std::endl;
	std::cout << "HostPort   = " << host_port << std::endl;
	std::cout << "LocalIp    = " << local_ip << std::endl;
	std::cout << "LocalPort  = " << local_port << std::endl;
	std::cout << "BufferIp   = " << buffer_ip << std::endl;
	std::cout << "BufferPort = " << buffer_port << std::endl << std::endl;

	// Define variables for the network socket descriptors
	//------------------------------------------------------------------------------
	socketd_t command_socket = 0;
	socketd_t listen_socket = 0;
	socketd_t data_socket = 0;

	try
	{
		if (buffer_ip == "-")
		{
			ft_server = start_server(buffer_port);
			std::cerr << "DEBUG: start_server  = " << ft_server << std::endl;
			if (ft_server!=0) exit(-1);
		}
		else
		{
			ft_server = open_connection(buffer_ip.c_str(), buffer_port);
			std::cerr << "DEBUG: open_connection  = " << ft_server << std::endl;
			if (ft_server<0) exit(-1);
		}

		int ret = InitNetworking();
		if (ret != NETWORKING_STARTUP_SUCCESS)
		return -1;

		socketd_t command_socket = 0;

		// Try to connect to the server; this port is used to send control
		// messages to the server
		//----------------------------------------------------------------------
		ret = EstablishNetworkConnection(&command_socket, host_ip, host_port);
		if (ret == SOCKET_ERROR)
		{
			CloseNetworkConnection(command_socket);
			CleanupNetworking();
			return -1;
		}

		// Setup the transceiver which handles basic communication and heartbeat
		// messages of the control pipe
		//----------------------------------------------------------------------
		Transceiver transceiver(command_socket);

		// Setup a message to retrieve the session id
		// (the session id is used as identifier in each message exchange)
		//----------------------------------------------------------------------
		std::string msg = SetupXMLMessage("0", CMD_GET_SESSION_ID);
		std::string reply = transceiver.SendReceive(msg);
		CheckServerReply(reply);
		std::string session_id = ParseXML(reply, GDS_XML_SESSION_ID_NODE);

		// Setup a message to retrieve the connected devices.
		// devices_in_use_only	'false'	derive information about all devices
		//						'true' derive information about devices
		//							   in use only
		//----------------------------------------------------------------------
		bool devices_in_use_only = false;
		std::string payload = EscapeXML(SetupXMLMessage(ToString(devices_in_use_only)));
		msg = SetupXMLMessage(session_id, CMD_GET_CONNECTED_DEVICES, payload);
		reply = transceiver.SendReceive(msg);
		CheckServerReply(reply);
		std::string device_info = EscapeXML(ParseXML(reply, GDS_XML_PAYLOAD_NODE), false);

		// Save the device list (for further usage)
		std::vector<Device> device_list = XMLToDevice(device_info);

		// Print the available devices
		if (device_list.empty())
		{
			std::cout << "No available devices on " + host_ip + "." << std::endl;
			return -1;
		}
		else
		{
			std::cout << "Available devices on " + host_ip + ":" << std::endl << std::endl;
			PrintDevices(&std::cout, device_list);
		}

		// Let the user choose a device
		//----------------------------------------------------------------------
		int device_index = 0;
		PromptUser(&std::cout, &std::cin, "Please select a device number: ", std::string(), &device_index, &CheckDeviceIndex, device_list.size() - 1);
		Device selected_device = device_list[device_index];

		// Let the user choose the duration of measurement and the sample rate
		//----------------------------------------------------------------------
		int duration = 0;
		int sample_rate = 0;
		PromptUser(&std::cout, &std::cin, "Please enter a valid sample rate (in Hz): ", std::string(), &sample_rate, &CheckSampleRate, selected_device);
		PromptUser(&std::cout, &std::cin, "Please enter the desired measurement time (in seconds): ", std::string(), &duration, &CheckPositive);

		std::cout << std::endl << "Device: " << selected_device.name_ << " @ session id: " << session_id << std::endl;
		std::cout << "Duration = " << duration << "s" << std::endl;
		std::cout << "Sample rate = " << sample_rate << "Hz" << std::endl;

		// Establish a socket which is listening for incoming connections
		//----------------------------------------------------------------------
		ListenOnNetwork(&listen_socket, local_port);
		if (listen_socket == SOCKET_ERROR)
		{
			CloseNetworkConnection(command_socket);
			CleanupNetworking();
			return -1;
		}

		// Tell the server where the client is accepting an incoming connection
		// (this connection is used to transfer measurement data)
		//----------------------------------------------------------------------
		std::pair<std::string, std::string> local_endpoint;
		local_endpoint.first = local_ip;
		local_endpoint.second = ToString(local_port);
		payload = EscapeXML(SetupXMLMessage(PairToXML(local_endpoint)));
		msg = SetupXMLMessage(session_id, CMD_SETUP_STREAMING, payload);
		reply = transceiver.SendReceive(msg);
		CheckServerReply(reply);

		// Accept the incoming connection (used to transfer measurement data)
		//----------------------------------------------------------------------
		data_socket = AcceptOnNetwork(listen_socket);
		if (data_socket == SOCKET_ERROR)
		{
			CloseNetworkConnection(listen_socket);
			CloseNetworkConnection(command_socket);
			CleanupNetworking();
			return -1;
		}

		// Open a data acquisition session
		//----------------------------------------------------------------------
		std::vector<std::string> device_name(1, selected_device.name_);
		payload = EscapeXML(SetupXMLMessage(VectorToXML(device_name)));
		msg = SetupXMLMessage(session_id, CMD_OPEN_DAQ_SESSION_EXCLUSIVELY, payload);
		reply = transceiver.SendReceive(msg);
		CheckServerReply(reply);

		// Setup the configuration (in XML) for the chosen device; the
		// configuration is stored in an extra variable (for further usage)
		//----------------------------------------------------------------------
		std::string config;
		switch (selected_device.type_)
		{
			case GHIAMP:
			config = SetupgHIampXMLConfig(selected_device.name_, ToString(sample_rate));
			break;
			case GUSBAMP:
			config = SetupgUSBampXMLConfig(selected_device.name_, ToString(sample_rate));
			break;
			case GNAUTILUS:
			config = SetupgNautilusXMLConfig(selected_device.name_, ToString(sample_rate));
			break;
			default:
			break;
		}

		payload = EscapeXML(config);
		msg = SetupXMLMessage(session_id, CMD_SET_CONFIGURATION, payload);
		reply = transceiver.SendReceive(msg);
		CheckServerReply(reply);

		// Retrieve information about the data (for one scan); the information
		// is stored in an extra variable (for further usage)
		//----------------------------------------------------------------------
		int scan_count = 1;
		payload = EscapeXML(SetupXMLMessage(ToString(scan_count)));
		msg = SetupXMLMessage(session_id, CMD_GET_DATA_INFO, payload);
		std::string data_info = transceiver.SendReceive(msg);

		// Start data acquisition and also wait for the event which indicates
		// that data acquisition has started
		//----------------------------------------------------------------------
		msg = SetupXMLMessage(session_id, CMD_START_ACQUISITION);
		reply = transceiver.SendReceive(msg, EVENT_DATA_ACQUISITION_STARTED);
		CheckServerReply(reply);

		// Start streaming data via network
		//----------------------------------------------------------------------
		msg = SetupXMLMessage(session_id, CMD_START_STREAMING);
		reply = transceiver.SendReceive(msg);
		CheckServerReply(reply);
		// Give the start streaming command some time (1s). This is a known
		// issue and will be changed in following releases.
		MKR_SLEEP(1);

		// Read measurement data; data will be saved to the file
		// data.bin
		//----------------------------------------------------------------------
		std::string sample_rate_node;
		switch (selected_device.type_)
		{
			case GHIAMP:
			sample_rate_node = XML_GHIAMP_SAMPLE_RATE_NODE;
			break;
			case GUSBAMP:
			sample_rate_node = XML_GUSBAMP_SAMPLE_RATE_NODE;
			break;
			case GNAUTILUS:
			sample_rate_node = XML_GNAUTILUS_SAMPLE_RATE_NODE;
			break;
			default:
			break;
		}

		ReadMeasurementData(&transceiver, data_socket, session_id, config, data_info, sample_rate_node, duration);

		// Stop streaming data via network
		//----------------------------------------------------------------------
		msg = SetupXMLMessage(session_id, CMD_STOP_STREAMING);
		reply = transceiver.SendReceive(msg);
		CheckServerReply(reply);

		// Stop data acquisition
		//----------------------------------------------------------------------
		msg = SetupXMLMessage(session_id, CMD_STOP_ACQUISITION);
		reply = transceiver.SendReceive(msg);
		CheckServerReply(reply);

		// Close data acquisition session
		//----------------------------------------------------------------------
		msg = SetupXMLMessage(session_id, CMD_CLOSE_DAQ_SESSION);
		reply = transceiver.SendReceive(msg);
		CheckServerReply(reply);

		// Disconnect from the server
		//----------------------------------------------------------------------
		msg = SetupXMLMessage(session_id, CMD_DISCONNECT);
		reply = transceiver.SendReceive(msg);
		CheckServerReply(reply);
	}
	catch (const int &error)
	{
		// Set the error code according to the retrieved error given by the
		// server's reply message
		ret_main = error;
	}

	// Close the socket by using the associated file descriptors and cleanup
	//------------------------------------------------------------------------------
	CloseNetworkConnection(command_socket);
	CloseNetworkConnection(data_socket);
	CloseNetworkConnection(listen_socket);
	CleanupNetworking();

	return ret_main;
}

std::string ToString(int v)
{
	std::stringstream converter;
	converter << v;
	return converter.str();
}

bool CheckDeviceIndex(int index, const size_t& max_index)
{
	if (index < 0)
	return false;

	return (size_t) index <= max_index;
}

bool CheckSampleRate(int sample_rate, const Device& device)
{
	return (0 < sample_rate);
}

bool CheckPositive(int number, const size_t& dummy)
{
	return (0 < number);
}

template<typename T1, typename T2>
void PromptUser(std::ostream* outstream, std::istream* instream, std::string message, std::string error_message, T1* result, bool (*constrain_function)(T1, const T2&), T2 constrain)
{
	std::string input;
	bool inputValid = false;

	do
	{
		*outstream << message;
		getline(*instream, input);

		// Convert from string to number safely.
		std::stringstream inputstream(input);
		inputValid = inputstream >> *result && (*constrain_function)(*result, constrain);

		if (!inputValid)
		{
			std::string error_string = (error_message.empty()) ? "Invalid value, please try again." : error_message;
			*outstream << error_string << std::endl;
		}
	} while (!inputValid);
}

void PrintDevices(std::ostream* outstream, const std::vector<Device>& device_list)
{
	*outstream << "    Nr.\tDevice name\tDevice type\tIn Use" << std::endl;
	*outstream << "    ------------------------------------------" << std::endl;

	for (std::vector<Device>::size_type i = 0; i != device_list.size(); i++)
	{
		*outstream << "    " << i;
		*outstream << "\t" << device_list[i].name_;
		*outstream << "\t" << device_list[i].full_type_;
		*outstream << "\t" << ((device_list[i].state_) ? "yes" : "no") << std::endl;
	}

	*outstream << std::endl << std::endl;
}
