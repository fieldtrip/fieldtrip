//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef PARSEXML_HPP_INCLUDED
#define PARSEXML_HPP_INCLUDED

#include <string>
#include <stdint.h>
#include <iostream>

#include "Device.hpp"
#include "EscapeXML.hpp"
#include "XMLDefinitions.hpp"
#include "XMLSTLConversions.hpp"

// This method returns the value which is stored within the XML element defined
// by the xml_node paramter.
//-------------------------------------------------------------------------------------
std::string ParseXML(std::string xml, std::string xml_node, size_t* pos_finally, size_t pos_begin);
std::string ParseXML(std::string xml, std::string xml_node);

// Retrieves the information about the devices connected to the server stored in an XML message.
//-------------------------------------------------------------------------------------
std::vector<struct Device> ParseXML(std::string xml, std::string xml_node_first, std::string xml_node_second, std::string xml_node_third);

// Retrieves the data stored in an XML message which stores the data info of
// measurement data
//-------------------------------------------------------------------------------------
void ParseXML(std::string xml, size_t* scan_count, size_t* channel_count, size_t* buffer_size_in_samples);

#endif // PARSEXML_HPP_INCLUDED
