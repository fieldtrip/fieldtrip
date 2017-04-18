//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef TRANSCEIVER_HPP_INCLUDED
#define TRANSCEIVER_HPP_INCLUDED

#include <ctime>
#include <string>
#include <list>
#include <iostream>

#include "gdsnetworkinterfacedemo_platform_specific.hpp"
#include "CheckServerReply.hpp"
#include "MetaHeader.hpp"
#include "Constants.hpp"
#include "ParseXML.hpp"

// Handles control message flow including heartbeat and acknowledge mechanism.
class Transceiver
{
public:
	// Initializes a new instance of the Transceiver class using the given socket and starts listening for incoming messages.
	explicit Transceiver(const socketd_t &socket);

	//Destructor. Stops receiving before resources are released.
	virtual ~Transceiver(void);

	// Sends the given command message and returns the received reply. The method blocks until a reply is received or the default timeout has elapsed.
	std::string SendReceive(const std::string& command_message, const std::string& event_code = "");

	// Sends the given command message only.
	bool Send(const std::string& command_message);

	// Waits for the given command to be received and returns true if so; returns false if the given timeout has elapsed or an error occurs. 
	// The method blocks until the given command is received or the given timeout has elapsed.
	bool Receive(const std::string& command_to_receive, const std::string& event_code_to_receive, std::string* received_reply, unsigned int timeout_milliseconds = 5000);

public:
	// Sends an acknowledge message for a received message represented by the given header.
	// Note: the class's member functions SendReceive, Send, and Receive manage acknowledgement automatically; they do not require the caller to use this method.
	static bool SendAcknowledge(const socketd_t &socket, const MetaHeader& received_message_header);

private:
	// Waits for an acknowledge message to be received for a sent message represented by the given header.
	bool WaitForAcknowledge(const MetaHeader& sent_message_header, unsigned int acknowledge_timeout_milliseconds = 5000);

	// Starts a thread that listens for incoming messages.
	void StartReceiving();

	// Stops listening for incoming messages.
	void StopReceiving();

private:
	// The thread handler that receives incoming messages in a loop until listening is stopped.
	static THREAD_CALLBACK_RETURN_TYPE THREAD_CALLBACK_CALLING_CONVENTION DoReceive(void *data);

private:
	const socketd_t &socket_;

	volatile bool receiving_;
	thread_t receiving_thread_;

	event_t message_received_;
	event_t acknowledge_received_;
	
	mutex_t received_messages_mutex_;
	mutex_t received_acknowledge_id_mutex_;

	std::list<std::string> received_messages_;
	uint64_t received_acknowledge_id_;
};

#endif // TRANSCEIVER_HPP_INCLUDED