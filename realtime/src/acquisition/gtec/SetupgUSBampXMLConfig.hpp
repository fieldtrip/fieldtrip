//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef SETUPGUSBAMPXMLCONFIG_HPP_INCLUDED
#define SETUPGUSBAMPXMLCONFIG_HPP_INCLUDED

#include <string>
#include "XMLDefinitions.hpp"
#include "XMLSTLConversions.hpp"

#define XML_GUSBAMP_PAIR_SECOND_CONFIG_META_INFO  " class_id=\"2\" tracking_level=\"0\" version=\"1\""
#define XML_GUSBAMP_BASE_META_INFO_BEGIN "<px class_id=\"4\" class_name=\"gAPI::gUSBamp::Config::Config\" tracking_level=\"1\" version=\"0\" object_id=\"_0\"><BaseConfig class_id=\"5\" tracking_level=\"0\" version=\"0\"><DeviceConfigurationBase class_id=\"3\" tracking_level=\"1\" version=\"0\" object_id=\"_1\"></DeviceConfigurationBase></BaseConfig>\n"
#define XML_GUSBAMP_BASE_META_INFO_END "</px>\n"
#define XML_GUSBAMP_VECTOR_DEFINITION_BEGIN "<TYPE_STD_VECTOR>\n"
#define XML_GUSBAMP_VECTOR_DEFINITION_END "</TYPE_STD_VECTOR>\n"

#define XML_GUSBAMP_BEGIN_CONFIG_SINGLE_DEVICE "<count>1</count>\n<item_version>0</item_version>\n<item class_id=\"1\" tracking_level=\"0\" version=\"0\">\n"
#define XML_GUSBAMP_END_CONFIG_SINGLE_DEVICE "</item>\n"

//___________________________________________________________________
// various defintions for the configuration
#define XML_GUSBAMP_SAMPLE_RATE_NODE "sample_rate_"
#define XML_GUSBAMP_SAMPLE_RATE_BEGIN "<sample_rate_>"
#define XML_GUSBAMP_SAMPLE_RATE_END "</sample_rate_>\n"
#define XML_GUSBAMP_NUMBER_OF_SCANS_BEGIN "<number_of_scans_>"
#define XML_GUSBAMP_NUMBER_OF_SCANS_END "</number_of_scans_>\n"
#define XML_GUSBAMP_SLAVE_MODE_BEGIN "<slave_>"
#define XML_GUSBAMP_SLAVE_MODE_END "</slave_>\n"
#define XML_GUSBAMP_SHORT_CUT_BEGIN "<short_cut_>"
#define XML_GUSBAMP_SHORT_CUT_END "</short_cut_>\n"
#define XML_GUSBAMP_DEVICE_MODE_BEGIN "<device_mode_>"
#define XML_GUSBAMP_DEVICE_MODE_END "</device_mode_>\n"
#define XML_GUSBAMP_SCAN_DIO_BEGIN "<scan_dio_>"
#define XML_GUSBAMP_SCAN_DIO_END "</scan_dio_>\n"

//___________________________________________________________________
// Begin analog out config definitions
#define XML_GUSBAMP_PARENT_ANALOG_OUT_BEGIN "<AnalogOut class_id=\"6\" tracking_level=\"0\" version=\"0\">\n"
#define XML_GUSBAMP_ANALOG_OUT_SHAPE_BEGIN "<shape_>"
#define XML_GUSBAMP_ANALOG_OUT_SHAPE_END "</shape_>\n"
#define XML_GUSBAMP_ANALOG_OUT_FREQUENCY_BEGIN "<frequency_>"
#define XML_GUSBAMP_ANALOG_OUT_FREQUENCY_END "</frequency_>\n"
#define XML_GUSBAMP_ANALOG_OUT_AMPLITUDE_BEGIN "<amplitude_>"
#define XML_GUSBAMP_ANALOG_OUT_AMPLITUDE_END "</amplitude_>\n"
#define XML_GUSBAMP_ANALOG_OUT_OFFSET_BEGIN "<offset_>"
#define XML_GUSBAMP_ANALOG_OUT_OFFSET_END "</offset_>\n"
#define XML_GUSBAMP_PARENT_ANALOG_OUT_END "</AnalogOut>\n"
//  End analog out config definitions
//___________________________________________________________________

//___________________________________________________________________
// definitions for common ground & reference
#define XML_GUSBAMP_PARENT_COMMON_GROUND_BEGIN "<common_ground_ class_id=\"7\" tracking_level=\"0\" version=\"0\">\n"
#define XML_GUSBAMP_PARENT_COMMON_GROUND_END "</common_ground_>\n"
#define XML_GUSBAMP_PARENT_COMMON_REFERENCE_BEGIN "<common_reference_>\n"
#define XML_GUSBAMP_PARENT_COMMON_REFERENCE_END "</common_reference_>\n"

//___________________________________________________________________
// definitions for the filter setttings
#define XML_GUSBAMP_PARENT_BANDPASS_BEGIN "<bandpass_>\n"
#define XML_GUSBAMP_PARENT_BANDPASS_END "</bandpass_>\n"
#define XML_GUSBAMP_PARENT_NOTCH_BEGIN "<notch_>\n"
#define XML_GUSBAMP_PARENT_NOTCH_END "</notch_>\n"

//___________________________________________________________________
// definition for the bipolar derivation
#define XML_GUSBAMP_PARENT_BIPOLAR_BEGIN "<bipolar_>\n"
#define XML_GUSBAMP_PARENT_BIPOLAR_END "</bipolar_>\n"

//___________________________________________________________________
// definitions for the measurement channels
#define XML_GUSBAMP_PARENT_ANALOG_IN_BEGIN "<analog_in_ class_id=\"11\" tracking_level=\"0\" version=\"0\">\n"
#define XML_GUSBAMP_PARENT_ANALOG_IN_END "</analog_in_>\n"

//___________________________________________________________________
// returns a g.USBamp config as XML string
std::string SetupgUSBampXMLConfig(std::string device_name,
	std::string sample_rate = "256",
	std::string number_of_scans = "0",
	std::string slave = "0",
	std::string short_cut = "0",
	std::string device_mode = "3", //1 = analog out, 3 = counter
	std::string scan_dio = "0",
	std::string analog_out_shape = "3",
	std::string analog_out_frequency = "10",
	std::string analog_out_amplitude = "100",
	std::string analog_out_offset = "0");

#endif // SETUPGUSBAMPXMLCONFIG
