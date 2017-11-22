//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef XMLDEFINITIONS_HPP_INCLUDED
#define XMLDEFINITIONS_HPP_INCLUDED

//____________________________________________________________________
// XML definitions from boost
#define GDS_XML_DECLARATION "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>\n"
#define GDS_XML_DOCTYPE "<!DOCTYPE boost_serialization>\n"
#define GDS_XML_METAINFO_BEGIN "<boost_serialization signature=\"serialization::archive\" version=\"9\">\n"
#define GDS_XML_METAINFO_END "</boost_serialization>\n"

#define GDS_XML_VALUE_NODE_BEGIN "<value class_id=\"0\" tracking_level=\"0\" version=\"0\">\n"
#define GDS_XML_VALUE_NODE_END "</value>\n"
#define GDS_XML_VALUE_NODE "value"

//___________________________________________________________________
// XML elemt definitions of a command message
#define GDS_XML_SESSION_ID_NODE "session_id_"
#define GDS_XML_COMMAND_NODE "command_"
#define GDS_XML_PAYLOAD_NODE "payload_"
#define GDS_XML_STATUS_NODE "status_"
#define GDS_XML_ERROR_NODE "error_"

//___________________________________________________________________
// XML elemt definitions of measurement data information
#define GDS_XML_DATA_INFO_NUMBER_OF_SCANS_NODE "number_of_scans_"
#define GDS_XML_DATA_INFO_CHANNELS_PER_DEVICE_NODE "channels_per_device_"
#define GDS_XML_DATA_INFO_BUFFER_SIZE_IN_SAMPLES_NODE "buffer_size_in_samples_"

//___________________________________________________________________
// XML elemt definitions of device information
#define GDS_XML_DEVICE_INFO_NAME "Name"
#define GDS_XML_DEVICE_INFO_TYPE "DeviceType"
#define GDS_XML_DEVICE_INFO_STATUS "InUse"

//____________________________________________________________________
// C++ STL to XML complients

// conversion to list (for example std::vector)
#define XML_STL_LIST_COUNT "count"
#define XML_STL_LIST_ITEM_VERSION "item_version"
#define XML_STL_LIST_ITEM "item"

// conversion to Device
#define XML_DEVICE_LIST_STATUS "InUse"
#define XML_DEVICE_LIST_NAME "Name"
#define XML_DEVICE_LIST_TYPE "DeviceType"

// conversion to std::pair
#define XML_STL_PAIR_FIRST "first"
#define XML_STL_PAIR_SECOND "second"

//____________________________________________________________________
// XML tokens and escape codes.
#define XML_TOKEN_QUOT "\""
#define XML_TOKEN_APOS "\'"
#define XML_TOKEN_LT "<"
#define XML_TOKEN_LTS "</"
#define XML_TOKEN_GT ">"
#define XML_TOKEN_AMP "&"

#define XML_ESCAPE_QUOT "&quot;"
#define XML_ESCAPE_APOS "&apos;"
#define XML_ESCAPE_LT "&lt;"
#define XML_ESCAPE_GT "&gt;"
#define XML_ESCAPE_AMP "&amp;"

#endif // XMLDEFINITIONS_HPP_INCLUDED
