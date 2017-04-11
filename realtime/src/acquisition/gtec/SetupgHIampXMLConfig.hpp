//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef SETUPGHIAMPXMLCONFIG_HPP_INCLUDED
#define SETUPGHIAMPXMLCONFIG_HPP_INCLUDED

#include <string>
#include "XMLDefinitions.hpp"
#include "XMLSTLConversions.hpp"

#define XML_GHIAMP_PAIR_SECOND_CONFIG_META_INFO  " class_id=\"2\" tracking_level=\"0\" version=\"1\""
#define XML_GHIAMP_BASE_META_INFO_BEGIN "<px class_id=\"4\" class_name=\"gAPI::gHIamp::Configuration\" tracking_level=\"1\" version=\"0\" object_id=\"_0\">\n<DeviceConfigurationBase class_id=\"3\" tracking_level=\"1\" version=\"0\" object_id=\"_1\"></DeviceConfigurationBase>\n<deviceConfiguration_ class_id=\"5\" tracking_level=\"0\" version=\"0\">\n"
#define XML_GHIAMP_BASE_META_INFO_END "</deviceConfiguration_></px>\n"
#define XML_GHIAMP_CHANNEL_VECTOR_ITEM_EXTENSION "<item class_id=\"7\" tracking_level=\"0\" version=\"0\">"
#define XML_GHIAMP_BEGIN_CONFIG_SINGLE_DEVICE "<count>1</count>\n<item_version>0</item_version>\n<item class_id=\"1\" tracking_level=\"0\" version=\"0\">\n"
#define XML_GHIAMP_END_CONFIG_SINGLE_DEVICE "</item>\n"

//___________________________________________________________________
// various defintions for the configuration
#define XML_GHIAMP_SAMPLE_RATE_NODE "SamplingRate"
#define XML_GHIAMP_SAMPLE_RATE_BEGIN "<SamplingRate>"
#define XML_GHIAMP_SAMPLE_RATE_END "</SamplingRate>\n"
#define XML_GHIAMP_NUMBER_OF_SCANS_BEGIN "<NumberOfScans>"
#define XML_GHIAMP_NUMBER_OF_SCANS_END "</NumberOfScans>\n"
#define XML_GHIAMP_TRIGGER_LINES_ENABLED_BEGIN "<TriggerLinesEnabled>"
#define XML_GHIAMP_TRIGGER_LINES_ENABLED_END "</TriggerLinesEnabled>\n"
#define XML_GHIAMP_HOLD_ENABLED_BEGIN "<HoldEnabled>"
#define XML_GHIAMP_HOLD_ENABLED_END "</HoldEnabled>\n"
#define XML_GHIAMP_SLAVE_MODE_BEGIN "<IsSlave>"
#define XML_GHIAMP_SLAVE_MODE_END "</IsSlave>\n"
#define XML_GHIAMP_COUNTER_ENABLED_BEGIN "<CounterEnabled>"
#define XML_GHIAMP_COUNTER_ENABLED_END "</CounterEnabled>\n"

#define XML_GHIAMP_PARENT_INTERNAL_SIGNAL_GENERATOR_BEGIN "<InternalSignalGenerator class_id=\"6\" tracking_level=\"0\" version=\"0\">\n"
#define XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_ENABLED_BEGIN "<Enabled>"
#define XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_ENABLED_END "</Enabled>\n"
#define XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_WAVE_SHAPE_BEGIN "<WaveShape>"
#define XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_WAVE_SHAPE_END "</WaveShape>\n"
#define XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_AMPLITUDE_BEGIN "<Amplitude>"
#define XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_AMPLITUDE_END "</Amplitude>"
#define XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_FREQUENCY_BEGIN "<Frequency>"
#define XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_FREQUENCY_END "</Frequency>\n"
#define XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_OFFSET_BEGIN "<Offset>"
#define XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_OFFSET_END "</Offset>\n"
#define XML_GHIAMP_PARENT_INTERANL_SIGNAL_GENERATOR_END "</InternalSignalGenerator>\n"

#define XML_GHIAMP_PARENT_CHANNELS_BEGIN "<Channels>\n"
#define XML_GHIAMP_CHANNEL_NUMBER_BEGIN "<ChannelNumber>"
#define XML_GHIAMP_CHANNEL_NUMBER_END "</ChannelNumber>\n"
#define XML_GHIAMP_CHANNEL_ACQUIRE_BEGIN "<Acquire>"
#define XML_GHIAMP_CHANNEL_ACQUIRE_END "</Acquire>\n"
#define XML_GHIAMP_CHANNEL_BANDPASS_FILTER_INDEX_BEGIN "<BandpassFilterIndex>"
#define XML_GHIAMP_CHANNEL_BANDPASS_FILTER_INDEX_END "</BandpassFilterIndex>\n"
#define XML_GHIAMP_CHANNEL_NOTCH_FILTER_INDEX_BEGIN "<NotchFilterIndex>"
#define XML_GHIAMP_CHANNEL_NOTCH_FILTER_INDEX_END "</NotchFilterIndex>\n"
#define XML_GHIAMP_CHANNEL_REFERENCE_BEGIN "<ReferenceChannel>"
#define XML_GHIAMP_CHANNEL_REFERENCE_END "</ReferenceChannel>\n"
#define XML_GHIAMP_PARENT_CHANNELS_END "</Channels>\n"

//___________________________________________________________________
// returns a g.USBamp config as XML string
std::string SetupgHIampXMLConfig(std::string device_name,
	std::string sampling_rate = "256", 
	size_t number_of_channels = 40,
	std::string number_of_scans = "0", 
	std::string trigger_lines_enabled = "0", 
	std::string hold_enabled = "0", 
	std::string is_slave = "0",
	std::string counter_enabled = "1",
	std::string signal_generator_enabled = "0",
	std::string wave_shape = "1",
	std::string amplitude = "7000",
	std::string frequency = "10",
	std::string offset = "-7000");

#endif // SETUPGHIAMPXMLCONFIG
