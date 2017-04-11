#include "Transceiver.hpp"


Transceiver::Transceiver(const socketd_t &socket)
	: socket_(socket),
	receiving_(false),
	receiving_thread_(0),
	message_received_(),
	acknowledge_received_(),
	received_messages_mutex_(),
	received_acknowledge_id_mutex_(),
	received_messages_(),
	received_acknowledge_id_(0)
	
{
	InitializeEvent(&message_received_);
	InitializeEvent(&acknowledge_received_);
	MKR_MUTEX_INIT(received_messages_mutex_);

	StartReceiving();
}


Transceiver::~Transceiver(void)
{
	StopReceiving();

	DestroyEvent(&message_received_);
	DestroyEvent(&acknowledge_received_);
	MKR_MUTEX_DESTROY(received_messages_mutex_);
}

std::string Transceiver::SendReceive(const std::string& command_message, const std::string& event_code)
{
	std::string reply;

    if (!Send(command_message))
        return "";

    if (!Receive(ParseXML(command_message, GDS_XML_COMMAND_NODE), event_code, &reply))
        return "";

    // if the received command matches the expectations then the server message is replied; otherwise, an empty string is replied
    return reply;
}

bool Transceiver::Send(const std::string& command_message)
{
	// retrieve the command from the given xml message
    std::string command(ParseXML(command_message, GDS_XML_COMMAND_NODE));

    if (command.empty())
        return false;

    // allocate memory for the buffer to be sent
    size_t buffer_size = command_message.size() + sizeof(MetaHeader);
    char* buffer = new char[buffer_size];

    // prepare the header
    MetaHeader header(command_message.size());

    // 1) copy the header to the beginning of the buffer
    // 2) append the payload to the buffer
    std::copy((char*) &header, (char*) &header + sizeof(MetaHeader), buffer);
    std::copy(command_message.begin(), command_message.end(), buffer + sizeof(MetaHeader));

    // write the whole header and payload to the socket
    int n = MKR_WRITE( socket_, buffer, buffer_size );

    delete[] buffer;
    buffer = 0;

    if (n == SOCKET_ERROR)
	{
        std::cerr << "ERROR: could not write to socketfd " << socket_ << std::endl;
        return false;
    }

    // check if the amount of data written is valid
    if (static_cast<size_t>(n) != buffer_size)
	{
        std::cerr << "ERROR: tried to write " << buffer_size << " bytes (MetaHeader + Payload)" << std::endl;
        return false;
    }

	// wait for acknowledge
	return WaitForAcknowledge(header);
}

bool Transceiver::Receive(const std::string& command_to_receive, const std::string& event_code_to_receive, std::string* received_reply, unsigned int timeout_milliseconds)
{
	bool command_match = false;
    bool event_match = false;
	std::string reply;
	std::list<std::string> pending_messages;

	time_t start_time_receive = time(0);

	do 
	{
		// wait for incoming message
		if (received_messages_.empty() && !WaitForEvent(&message_received_, timeout_milliseconds))
		{
			std::cerr << "ERROR: waiting for reply from server timed out, no reply received" << std::endl;
			return false;
		}
		 
		// make a local copy of the message list to avoid race conditions
		MKR_MUTEX_LOCK(received_messages_mutex_);
		pending_messages = received_messages_;
		MKR_MUTEX_UNLOCK(received_messages_mutex_);

		// check if one of the received messages is the dedicated message
		for (std::list<std::string>::iterator it = pending_messages.begin(); it != pending_messages.end(); ++it)
		{				
			std::string server_response = *it;

			// retrieve the received command; if the command does not match the expected command ( = the command sent to the server before ) then
			// read again from the socket (the server may have emitted asynchronous event messages which are received here as well)
			std::string command_replied(ParseXML(server_response, GDS_XML_COMMAND_NODE));

			if (command_replied == EVENT_DATA_ACQUISITION_ERROR)
			{
				CheckServerReply(server_response, true, false);
				continue;
			}

			// loop until the reply for the given command (and if selected, the given event as well) has been received
			if (!command_match)
				command_match = (command_replied == command_to_receive);

			if (!event_match)
				event_match = (command_replied == event_code_to_receive);

			if (command_match)
			{
				if (reply.empty())
					reply = server_response;

				// if the received command matches the expectations (and the event as well if specified) we've achieved our goal
				if (event_code_to_receive.empty() || event_match)
				{
					MKR_MUTEX_LOCK(received_messages_mutex_);
					received_messages_.remove(*it);
					MKR_MUTEX_UNLOCK(received_messages_mutex_);

					*received_reply = reply;
					return true;
				}
			}
		}
	} while (difftime(time(0), start_time_receive) < (timeout_milliseconds / 1000.0));

    std::cerr << "ERROR: waiting for reply from server timed out, no reply received" << std::endl;
    return false;
}

bool Transceiver::SendAcknowledge(const socketd_t &socket, const MetaHeader& received_message_header)
{
	MetaHeader acknowledge_header;
    acknowledge_header.acknowledge_ = 1;
    acknowledge_header.package_id_ = received_message_header.package_id_;

    // confirm data reception with acknowledge message (header only with acknowledge flag set)
    int n = MKR_WRITE(socket, (char*) &acknowledge_header, sizeof(MetaHeader));

    if (n == SOCKET_ERROR)
	{
        std::cerr << "ERROR: could not write acknowledge to socket " << socket << std::endl;
        return false;
    }

    if (static_cast<size_t>(n) != sizeof(MetaHeader))
	{
        std::cerr << "ERROR: tried to write " << sizeof(MetaHeader) << " bytes acknowledge (MetaHeader)" << std::endl;
        return false;
    }

    return true;
}

bool Transceiver::WaitForAcknowledge(const MetaHeader& sent_message_header, unsigned int acknowledge_timeout_milliseconds)
{
	// wait for acknowledge
	if (!WaitForEvent(&acknowledge_received_, 5000))
	{
		std::cerr << "ERROR: waiting for acknowledge timed out, no acknowledge received" << std::endl;
		return false;
	}

	MKR_MUTEX_LOCK(received_acknowledge_id_mutex_);
	uint64_t received_acknowledge_id = received_acknowledge_id_;
	MKR_MUTEX_UNLOCK(received_acknowledge_id_mutex_);

	if (received_acknowledge_id != sent_message_header.package_id_)
	{
		std::cerr << "ERROR: acknowledge received for different message" << std::endl;
		return false;
	}
       
	return true;
}

void Transceiver::StartReceiving()
{
	//check if thread is not running already
	if (receiving_ || receiving_thread_ != 0)
		return;

	receiving_ = MKR_THREAD_CREATE(&receiving_thread_, &DoReceive, (void*) this);
}

void Transceiver::StopReceiving()
{
	receiving_ = false;

	// if a disconnect request has been sent to the server beforehand, the server 
	// closes the connection gracefully which forces MKR_READ to return with zero
	MKR_THREAD_JOIN(receiving_thread_);
	MKR_THREAD_DESTROY(receiving_thread_);
	receiving_thread_ = 0;
}

THREAD_CALLBACK_RETURN_TYPE THREAD_CALLBACK_CALLING_CONVENTION Transceiver::DoReceive(void *data)
{
	Transceiver *transceiver = reinterpret_cast<Transceiver*>(data);
	MetaHeader received_header;

	while (transceiver->receiving_)
	{
		// read header (every GDS message starts with a header)
		int n = MKR_READ(transceiver->socket_, (char*) &received_header, sizeof(MetaHeader));

		// check if shutdown is pending; if a disconnect request has been sent to the server beforehand,
		// the server closes the connection gracefully which forces MKR_READ to return with zero
		if (!transceiver->receiving_)
			break;

		if (n == SOCKET_ERROR)
		{
			std::cerr << "ERROR: could not read header from socket " << transceiver->socket_ << std::endl;
			continue;
		}

		// check if the size of read data is valid
		if (static_cast<size_t>(n) != sizeof(MetaHeader))
		{
			std::cerr << "ERROR: tried to read " << sizeof(MetaHeader) << " bytes (MetaHeader), read n = " << n << " bytes" << std::endl;
			continue;
		}

		// check if the header identifier is valid
		if (received_header.header_identifier_ != GDSMSG)
		{
			std::cerr << "ERROR: header identifier is not valid" << std::endl;
			continue;
		}

		// receive payload if any (heartbeat and acknowledge come without payload)
		if (received_header.size_ > 0)
		{
			try
			{
				// prepare the buffer to save the payload
				char *buffer = new char[received_header.size_];

				// read the expected amount of data ( given by the header )
				int n = MKR_READ(transceiver->socket_, buffer, received_header.size_);

				while (static_cast<unsigned int>(n) != received_header.size_)
					n += MKR_READ(transceiver->socket_, buffer + n, received_header.size_ - n);

				//
				// TODO: To differ between ordinary replies and asynchronous event messages from the server,
				// here would be a good place to validate the correctness of the message, sort it out, and
				// raise a (non-blocking) callback if an asynchronous event message has been received.
				// Ordinary replies should be processed as follows.
				//

				// store received message
				MKR_MUTEX_LOCK(transceiver->received_messages_mutex_);
				transceiver->received_messages_.push_back(std::string(buffer, received_header.size_));
				MKR_MUTEX_UNLOCK(transceiver->received_messages_mutex_);

				// notify pending user that a message has been received
				SetEvent(&transceiver->message_received_);

				delete[] buffer;
			} 
			catch (const std::bad_alloc&) 
			{
				std::cerr << "ERROR: could not allocate memory of size received_header.size_ = " << std::hex << received_header.size_ << std::endl;
				continue;
			}
		}

		// acknowledge a common or heartbeat message that has not set acknowledge flag
		if (received_header.acknowledge_ == 0)
			Transceiver::SendAcknowledge(transceiver->socket_, received_header);
		else
		{
			// if an acknowledge message has been received, notify the send method about it
			MKR_MUTEX_LOCK(transceiver->received_acknowledge_id_mutex_);
			transceiver->received_acknowledge_id_ = received_header.package_id_;
			MKR_MUTEX_UNLOCK(transceiver->received_acknowledge_id_mutex_);

			SetEvent(&transceiver->acknowledge_received_);
		}
	}

	MKR_THREAD_EXIT();
	return 0;
}