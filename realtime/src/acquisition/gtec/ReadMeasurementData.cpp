//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#include "ReadMeasurementData.hpp"
#include "interface.h"

// FieldTrip buffer
extern int ft_server;

void ReadMeasurementData(Transceiver *control_command_transceiver,
	socketd_t data_socketfd,
	std::string session_id,
	std::string xml_config,
	std::string xml_data_info,
	std::string sample_rate_parent_node,
	unsigned int duration,
	std::string data_file,
	bool disp)
	{
		FILE *file = fopen(data_file.c_str(), "wb");

		if (file == 0)
		{
			std::cerr << "ERROR: could not create or open file '" << data_file << "' for writing." << std::endl;
			return;
		}

		size_t scan_count = 0;
		size_t channels_count = 0;
		size_t buffer_size_per_scan = 0;
		ParseXML(xml_data_info, &scan_count, &channels_count, &buffer_size_per_scan);

		//std::cerr << "DEBUG: scan_count             = " << scan_count << std::endl;
		//std::cerr << "DEBUG: channels_count         = " << channels_count << std::endl;
		//std::cerr << "DEBUG: buffer_size_per_scan   = " << buffer_size_per_scan << std::endl;

		scan_count = 0;
		uint64_t buffer_size_seconds = 5;
		uint64_t sample_rate = atoi(ParseXML(xml_config, sample_rate_parent_node).c_str());
		size_t buffer_size_in_samples = buffer_size_seconds * sample_rate * buffer_size_per_scan;
		float* buffer = new float[buffer_size_in_samples];

		// std::cerr << "DEBUG: buffer_size_seconds    = " << buffer_size_seconds << std::endl;
		// std::cerr << "DEBUG: sample_rate            = " << sample_rate << std::endl;
		// std::cerr << "DEBUG: buffer_size_in_samples = " << buffer_size_in_samples << std::endl;

		int ft_status = write_header(ft_server, DATATYPE_FLOAT32, channels_count, sample_rate);
		// std::cerr << "DEBUG: write_header = " << ft_status << std::endl;
		if (ft_status!=0) exit(-1);

		uint64_t total_acquired_scans = 0;
		uint64_t total_scans_to_acquire = duration * sample_rate;

		if (disp)
		std::cout << std::endl << "Start reading measurement data: Expect about " << total_scans_to_acquire << " scans" << std::endl;

		while (total_acquired_scans < total_scans_to_acquire)
		{
			std::string data_info = SetupXMLMessage(scan_count, 0, buffer_size_in_samples);
			std::string cmd = SetupXMLMessage(session_id, CMD_GET_DATA, EscapeXML(data_info));
			std::string reply;

			// send data request
			if (!control_command_transceiver->Send(cmd))
			continue;

			// receive data header
			MetaHeader data_header;
			size_t bytes_read = MKR_READ( data_socketfd, (char*) &data_header, sizeof( MetaHeader ) );

			if (bytes_read != sizeof(MetaHeader))
			{
				std::cerr << "ERROR: was not able to read the header. read " << bytes_read << " instead of  " << sizeof(MetaHeader) << std::endl;
				continue;
			}

			// receive data payload
			uint64_t transfer_parts_count = (uint64_t) floor(data_header.size_ / (double) MAX_TRANSFER_SIZE);
			uint64_t last_transfer_part_size = data_header.size_ % MAX_TRANSFER_SIZE;

			for (uint64_t i = 0; i < transfer_parts_count + 1; i++)
			{
				uint64_t transfer_part_size = (i < transfer_parts_count) ? MAX_TRANSFER_SIZE : last_transfer_part_size;
				bytes_read = 0;

				while (bytes_read != transfer_part_size)
				{
					int ret = MKR_READ( data_socketfd, (char*) buffer + i*MAX_TRANSFER_SIZE + bytes_read, transfer_part_size - bytes_read );
					if (ret == SOCKET_ERROR)
					std::cout << "ERROR: on reading data from the socket" << std::endl;

					bytes_read += ((ret < 0) ? 0 : ret);
				}
			}

			// confirm data reception with acknowledge message (header only with acknowledge flag set)
			if (!Transceiver::SendAcknowledge(data_socketfd, data_header))
			continue;

			// receive reply for data request
			if (!control_command_transceiver->Receive(CMD_GET_DATA, "", &reply))
			reply = "";

			try
			{
				CheckServerReply(reply);
			}
			catch (...)
			{
				delete[] buffer;
				buffer = 0;
				throw;
			}

			std::string payload = EscapeXML(ParseXML(reply, GDS_XML_PAYLOAD_NODE), false);
			size_t scans_available = atoi(ParseXML(payload, GDS_XML_VALUE_NODE).c_str());

			total_acquired_scans += scans_available;

			if (scans_available > 0) {
				// std::cerr << "DEBUG: scans_available = " << scans_available << std::endl;
				fwrite((void*) buffer, sizeof(float), scans_available * buffer_size_per_scan, file);
				ft_status = write_data(ft_server, DATATYPE_FLOAT32, channels_count, scans_available * buffer_size_per_scan / channels_count, (void *)buffer);
				// std::cerr << "DEBUG: write_data  = " << ft_status << std::endl;
				if (ft_status!=0) exit(-1);
			}

			if (disp)
			std::cout << TERMINAL_CARRIAGE_RETURN_ESCAPE_CODE << total_acquired_scans << " / " << total_scans_to_acquire << " scans acquired" << std::flush;
		}

		ft_status = close_connection(ft_server);
		// std::cerr << "DEBUG: close_connection  = " << ft_status << std::endl;
		if (ft_status!=0) exit(-1);

		if (disp)
		std::cout << std::endl;

		delete[] buffer;
		buffer = 0;

		if (file != 0)
		{
			fclose(file);
			file = 0;
		}
	}
