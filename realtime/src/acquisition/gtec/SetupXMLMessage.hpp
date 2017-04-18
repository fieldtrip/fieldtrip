//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef SETUPXMLMESSAGE_HPP_INCLUDED
#define SETUPXMLMESSAGE_HPP_INCLUDED

#include <string>
#include <sstream>

#include "Constants.hpp"
#include "XMLDefinitions.hpp"
#include "XMLSTLConversions.hpp"

std::string SetupXMLElement(std::string xml_node, std::string value, bool new_line = false);
std::string SetupXMLElement(std::string xml_node, size_t value, bool new_line = false);
std::string SetupXMLMessage(std::string session_id, std::string command, std::string payload = "");
std::string SetupXMLMessage(std::string xml_payload);
std::string SetupXMLMessage(size_t scan_count, size_t channel_count, size_t buffer_size_in_samples);

#endif // SETUPXMLMESSAGE_HPP_INCLUDED
