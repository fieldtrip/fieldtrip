//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef SETUPGNAUTILUSXMLCONFIG_HPP_INCLUDED
#define SETUPGNAUTILUSXMLCONFIG_HPP_INCLUDED

#include <string>
#include "XMLDefinitions.hpp"
#include "XMLSTLConversions.hpp"

#define XML_GNAUTILUS_PAIR_SECOND_CONFIG_META_INFO " class_id=\"2\" tracking_level=\"0\" version=\"1\""
#define XML_GNAUTILUS_BASE_META_INFO_BEGIN "<px class_id=\"4\" class_name=\"gAPI::gNautilus::Configuration\" tracking_level=\"1\" version=\"0\" object_id=\"_0\"><DeviceConfigurationBase class_id=\"3\" tracking_level=\"1\" version=\"0\" object_id=\"_1\"></DeviceConfigurationBase>"
#define XML_GNAUTILUS_BASE_META_INFO_END "</px>\n"
#define XML_GNAUTILUS_CHANNEL_VECTOR_ITEM_EXTENSION "<item_version>0</item_version>\n<item class_id=\"6\" tracking_level=\"0\" version=\"0\">\n"
#define XML_GNAUTILUS_BEGIN_CONFIG_SINGLE_DEVICE "<count>1</count><item_version>0</item_version><item class_id=\"1\" tracking_level=\"0\" version=\"0\">\n"
#define XML_GNAUTILUS_END_CONFIG_SINGLE_DEVICE "</item>\n"

//___________________________________________________________________
// various defintions for the configuration
#define XML_GNAUTILUS_SAMPLE_RATE_NODE "sample_rate_"
#define XML_GNAUTILUS_SAMPLE_RATE_BEGIN "<sample_rate_>"
#define XML_GNAUTILUS_SAMPLE_RATE_END "</sample_rate_>\n"
#define XML_GNAUTILUS_SLAVE_MODE_BEGIN "<slave_>"
#define XML_GNAUTILUS_SLAVE_MODE_END "</slave_>\n"
#define XML_GNAUTILUS_INPUT_SIGNAL_BEGIN "<input_signal_>"
#define XML_GNAUTILUS_INPUT_SIGNAL_END "</input_signal_>\n"
#define XML_GNAUTILUS_NUMBER_OF_SCANS_BEGIN "<number_of_scans_>"
#define XML_GNAUTILUS_NUMBER_OF_SCANS_END "</number_of_scans_>\n"
#define XML_GNAUTILUS_NOISE_REDUCTION_BEGIN "<noise_reduction_>"
#define XML_GNAUTILUS_NOISE_REDUCTION_END "</noise_reduction_>\n"
#define XML_GNAUTILUS_CAR_BEGIN "<car_>"
#define XML_GNAUTILUS_CAR_END "</car_>\n"
#define XML_GNAUTILUS_ACCELERATION_DATA_BEGIN "<acceleration_data_>"
#define XML_GNAUTILUS_ACCELERATION_DATA_END "</acceleration_data_>\n"
#define XML_GNAUTILUS_COUNTER_BEGIN "<counter_>"
#define XML_GNAUTILUS_COUNTER_END "</counter_>\n"
#define XML_GNAUTILUS_LINK_QUALITY_INFORMATION_BEGIN "<link_quality_information_>"
#define XML_GNAUTILUS_LINK_QUALITY_INFORMATION_END "</link_quality_information_>\n"
#define XML_GNAUTILUS_BATTERY_LEVEL_BEGIN "<battery_level_>"
#define XML_GNAUTILUS_BATTERY_LEVEL_END "</battery_level_>\n"
#define XML_GNAUTILUS_DIGITAL_IOS_BEGIN "<digital_ios_>"
#define XML_GNAUTILUS_DIGITAL_IOS_END "</digital_ios_\n>"
#define XML_GNAUTILUS_VALIDATION_INDICATOR_BEGIN "<validation_indicator_>"
#define XML_GNAUTILUS_VALIDATION_INDICATOR_END "</validation_indicator_>\n"
#define XML_GNAUTILUS_NETWORK_CHANNEL_BEGIN "<network_channel_>"
#define XML_GNAUTILUS_NETWORK_CHANNEL_END "</network_channel_>\n"

#define XML_GNAUTILUS_PARENT_CHANNELS_BEGIN "<channel_configurations_ class_id=\"5\" tracking_level=\"0\" version=\"0\">\n"
#define XML_GNAUTILUS_CHANNEL_NUMBER_BEGIN  "<channel_number_>"
#define XML_GNAUTILUS_CHANNEL_NUMBER_END "</channel_number_>\n"
#define XML_GNAUTILUS_CHANNEL_ENABLED_BEGIN "<enabled_>"
#define XML_GNAUTILUS_CHANNEL_ENABLED_END "</enabled_>\n"
#define XML_GNAUTILUS_CHANNEL_SENSITIVITY_BEGIN "<sensitivity_>"
#define XML_GNAUTILUS_CHANNEL_SENSITIVITY_END "</sensitivity_>\n"
#define XML_GNAUTILUS_CHANNEL_USED_FOR_NOISE_REDUCTION_BEGIN "<used_for_noise_reduction_>"
#define XML_GNAUTILUS_CHANNEL_USED_FOR_NOISE_REDUCTION_END "</used_for_noise_reduction_>\n"
#define XML_GNAUTILUS_CHANNEL_USED_FOR_CAR_BEGIN "<used_for_car_>"
#define XML_GNAUTILUS_CHANNEL_USED_FOR_CAR_END "</used_for_car_>\n"
#define XML_GNAUTILUS_CHANNEL_BANDPASS_FILTER_INDEX_BEGIN "<bandpass_filter_index_>"
#define xML_GNAUTILUS_CHANNEL_BANDPASS_FILTER_INDEX_END "</bandpass_filter_index_>\n"
#define XML_GNAUTILUS_CHANNEL_NOTCH_FILTER_INDEX_BEGIN "<notch_filter_index_>"
#define XML_GNAUTILUS_CHANNEL_NOTCH_FILTER_INDEX_END "</notch_filter_index_>\n"
#define XML_GNAUTILUS_CHANNEL_BIPOLAR_CHANNEL_BEGIN "<bipolar_channel_>"
#define XML_GNAUTILUS_CHANNEL_BIPOLAR_CHANNEL_END "</bipolar_channel_>\n"
#define XML_GNAUTILUS_PARENT_CHANNELS_END "</channel_configurations_>\n"

//___________________________________________________________________
// returns a g.Nautilus config as XML string
std::string SetupgNautilusXMLConfig(std::string device_name,
	std::string sampling_rate = "250",
	size_t number_of_channels = 32,
	std::string slave = "0",
	std::string number_of_scans = "0",
	std::string input_signal = "0",
	std::string acceleration_data = "0",
	std::string counter = "0",
	std::string link_quality_information = "0",
	std::string battery_level = "0",
	std::string digital_ios = "0",
	std::string validation_indicator = "0",
	std::string network_channel = "19");

#endif // SETUPGNAUTILUSXMLCONFIG
