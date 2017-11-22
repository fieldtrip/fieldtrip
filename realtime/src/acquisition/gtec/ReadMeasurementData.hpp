//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef READMEASUREMENTDATA_HPP_INCLUDED
#define READMEASUREMENTDATA_HPP_INCLUDED

#include <math.h>
#include <cstdio>

#include "SetupgUSBampXMLConfig.hpp"
#include "ParseXML.hpp"
#include "SetupXMLMessage.hpp"
#include "Transceiver.hpp"
#include "CheckServerReply.hpp"

#define TERMINAL_CLEAR_TO_THE_LEFT_ESCAPE_CODE "\0"
#define TERMINAL_CARRIAGE_RETURN_ESCAPE_CODE "\r"
#define MAX_TRANSFER_SIZE (64*1024 - 1)

// This method is responsible for the whole process of reading the
// measurement data transfered via the network
//-------------------------------------------------------------------------------------
void ReadMeasurementData(Transceiver *control_command_transceiver,
	socketd_t data_socketfd,
	std::string session_id, 
	std::string xml_config,
	std::string xml_data_info,
	std::string sample_rateparent_node,
	unsigned int duration = 10,
	std::string data_file = std::string("data.bin"), 
	bool disp = true);

#endif // READMEASUREMENTDATA_HPP_INCLUDED
