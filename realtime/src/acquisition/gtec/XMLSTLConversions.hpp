//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef XMLSTLCONVERSIONS_HPP_INCLUDED
#define XMLSTLCONVERSIONS_HPP_INCLUDED

#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>

#include "XMLDefinitions.hpp"
#include "ParseXML.hpp"
#include "Device.hpp"

static const std::string xml_stl_list_count(XML_STL_LIST_COUNT);
static const std::string token_count_begin("<" + xml_stl_list_count + ">");
static const size_t size_token_count_begin(token_count_begin.size());
static const std::string token_count_end("</" + xml_stl_list_count + ">\n");
static const size_t size_token_count_end(token_count_end.size());
static const std::string xml_stl_list_item_version(XML_STL_LIST_ITEM_VERSION);
static const std::string xml_stl_list_item(XML_STL_LIST_ITEM);
static const std::string token_item_begin("<" + xml_stl_list_item + ">");
static const size_t size_token_begin(token_item_begin.size());
static const std::string token_item_end("</" + xml_stl_list_item + ">\n");
static size_t size_token_end(token_item_end.size());

static const std::string xml_stl_pair_first(XML_STL_PAIR_FIRST);
static const std::string token_first_begin("<" + xml_stl_pair_first + ">");
static const size_t size_token_first_begin(token_first_begin.size());
static const std::string token_first_end("</" + xml_stl_pair_first + ">\n");
static const std::string xml_stl_pair_second(XML_STL_PAIR_SECOND);
static const std::string token_second_begin("<" + xml_stl_pair_second + ">");
static const size_t size_token_second_begin(token_second_begin.size());
static const std::string token_second_end("</" + xml_stl_pair_second + ">\n");

static const std::string xml_device_list_status(XML_DEVICE_LIST_STATUS);
static const std::string token_status_begin("<" + xml_device_list_status + ">");
static const size_t size_token_status_begin(token_count_begin.size());
static const std::string token_status_end("</" + xml_device_list_status + ">\n");
static const size_t size_token_status_end(token_status_end.size());
static const std::string xml_device_list_name(XML_DEVICE_LIST_NAME);
static const std::string token_name_begin("<" + xml_device_list_name + ">");
static const size_t size_token_name_begin(token_name_begin.size());
static const std::string token_name_end("</" + xml_device_list_name + ">\n");
static const size_t size_token_name_end(token_name_end.size());
static const std::string xml_device_list_type(XML_DEVICE_LIST_TYPE);
static const std::string token_type_begin("<" + xml_device_list_type + ">");
static const size_t size_token_type_begin(token_type_begin.size());
static const std::string token_type_end("</" + xml_device_list_type + ">\n");
static const size_t size_token_type_end(token_type_end.size());

// Builds a std vector from the given xml std vector representation.
// pos_finally inidicates where the last item has been found within xml.
// pos_init tells where to start.
//------------------------------------------------------------------------------
std::vector<std::string> XMLToVector(std::string xml, size_t* pos_finally, size_t pos_init = 0);

// Builds a std vector from the given xml std vector representation.
//------------------------------------------------------------------------------
std::vector<std::string> XMLToVector(std::string xml);

// Builds a std vector of Devices from the given xml std vector representation.
//------------------------------------------------------------------------------
std::vector<struct Device> XMLToDevice(std::string xml);

// Returns a XML string which represents the given vector. Additonally, the item 
// version field can be disabled.
//------------------------------------------------------------------------------
std::string VectorToXML(std::vector<std::string> v, bool enable_item_version = true, std::string item_version = "0");

// Returns a XML string which represents the given vector. Additonally, the item
// can be specified.
//------------------------------------------------------------------------------
std::string VectorToXML(std::vector<std::string> v, std::string item_extension);

// Returns a std pair represented by the XML string.
//------------------------------------------------------------------------------
std::pair<std::string, std::string> XMLToPair(std::string xml);

// Returns a XML string which represents the given pair.
//------------------------------------------------------------------------------
std::string PairToXML(std::pair<std::string, std::string> p, std::string xml_meta_info_first = "", std::string xml_meta_info = "");

#endif // XMLVECTORCONVERSIONS_HPP_INCLUDED
