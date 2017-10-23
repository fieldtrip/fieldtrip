#include "SetupgNautilusXMLConfig.hpp"

std::string SetupgNautilusXMLConfig(std::string device_name,
        std::string sampling_rate, 
		size_t number_of_channels, 
		std::string slave,
        std::string number_of_scans, 
		std::string input_signal,
        std::string acceleration_data, 
		std::string counter,
        std::string link_quality_information, 
		std::string battery_level,
        std::string digital_ios, 
		std::string validation_indicator,
        std::string network_channel) 
{
    std::string xml;

    xml = GDS_XML_DECLARATION;
    xml += GDS_XML_DOCTYPE;
    xml += GDS_XML_METAINFO_BEGIN;
    xml += GDS_XML_VALUE_NODE_BEGIN;
    xml += XML_GNAUTILUS_BEGIN_CONFIG_SINGLE_DEVICE;

    std::string gnautilus_config(XML_GNAUTILUS_BASE_META_INFO_BEGIN);

    gnautilus_config += XML_GNAUTILUS_SLAVE_MODE_BEGIN + slave + XML_GNAUTILUS_SLAVE_MODE_END;
    gnautilus_config += XML_GNAUTILUS_INPUT_SIGNAL_BEGIN + input_signal + XML_GNAUTILUS_INPUT_SIGNAL_END;
    gnautilus_config += XML_GNAUTILUS_NUMBER_OF_SCANS_BEGIN + number_of_scans + XML_GNAUTILUS_NUMBER_OF_SCANS_END;
    gnautilus_config += XML_GNAUTILUS_SAMPLE_RATE_BEGIN + sampling_rate + XML_GNAUTILUS_SAMPLE_RATE_END;
    gnautilus_config += XML_GNAUTILUS_NOISE_REDUCTION_BEGIN + std::string("0") + XML_GNAUTILUS_NOISE_REDUCTION_END;
    gnautilus_config += XML_GNAUTILUS_CAR_BEGIN + std::string("0") + XML_GNAUTILUS_CAR_END;
    gnautilus_config += XML_GNAUTILUS_ACCELERATION_DATA_BEGIN + acceleration_data + XML_GNAUTILUS_ACCELERATION_DATA_END;
    gnautilus_config += XML_GNAUTILUS_COUNTER_BEGIN + counter + XML_GNAUTILUS_COUNTER_END;
    gnautilus_config += XML_GNAUTILUS_LINK_QUALITY_INFORMATION_BEGIN + link_quality_information + XML_GNAUTILUS_LINK_QUALITY_INFORMATION_END;
    gnautilus_config += XML_GNAUTILUS_BATTERY_LEVEL_BEGIN + battery_level + XML_GNAUTILUS_BATTERY_LEVEL_END;
    gnautilus_config += XML_GNAUTILUS_DIGITAL_IOS_BEGIN + digital_ios + XML_GNAUTILUS_DIGITAL_IOS_END;
    gnautilus_config += XML_GNAUTILUS_VALIDATION_INDICATOR_BEGIN + validation_indicator + XML_GNAUTILUS_VALIDATION_INDICATOR_END;
    gnautilus_config += XML_GNAUTILUS_NETWORK_CHANNEL_BEGIN + network_channel + XML_GNAUTILUS_NETWORK_CHANNEL_END;

    gnautilus_config += XML_GNAUTILUS_PARENT_CHANNELS_BEGIN;

    std::vector<std::string> channels;
    for (size_t i = 0; i < 64; i++) 
	{
        std::stringstream sstr;
        sstr << i;
        std::string channel_setup;
        channel_setup = XML_GNAUTILUS_CHANNEL_NUMBER_BEGIN + sstr.str() + XML_GNAUTILUS_CHANNEL_NUMBER_END;
        sstr.str("");
        sstr << (i < number_of_channels);
        channel_setup += XML_GNAUTILUS_CHANNEL_ENABLED_BEGIN + sstr.str() + XML_GNAUTILUS_CHANNEL_ENABLED_END;
        channel_setup += XML_GNAUTILUS_CHANNEL_SENSITIVITY_BEGIN + std::string("187500") + XML_GNAUTILUS_CHANNEL_SENSITIVITY_END;
        channel_setup += XML_GNAUTILUS_CHANNEL_USED_FOR_NOISE_REDUCTION_BEGIN + std::string("0") + XML_GNAUTILUS_CHANNEL_USED_FOR_NOISE_REDUCTION_END;
        channel_setup += XML_GNAUTILUS_CHANNEL_USED_FOR_CAR_BEGIN + std::string("1") + XML_GNAUTILUS_CHANNEL_USED_FOR_CAR_END;
        channel_setup += XML_GNAUTILUS_CHANNEL_BANDPASS_FILTER_INDEX_BEGIN + std::string("-1") + xML_GNAUTILUS_CHANNEL_BANDPASS_FILTER_INDEX_END;
        channel_setup += XML_GNAUTILUS_CHANNEL_NOTCH_FILTER_INDEX_BEGIN + std::string("-1") + XML_GNAUTILUS_CHANNEL_NOTCH_FILTER_INDEX_END;
        channel_setup += XML_GNAUTILUS_CHANNEL_BIPOLAR_CHANNEL_BEGIN + std::string("-1") + XML_GNAUTILUS_CHANNEL_BIPOLAR_CHANNEL_END;
        channels.push_back(channel_setup);
    }

    gnautilus_config += VectorToXML(channels, std::string(XML_GNAUTILUS_CHANNEL_VECTOR_ITEM_EXTENSION));
    gnautilus_config += XML_GNAUTILUS_PARENT_CHANNELS_END;

    gnautilus_config += XML_GNAUTILUS_BASE_META_INFO_END;

    std::pair<std::string, std::string> device2config;
    device2config.first = device_name;
    device2config.second = gnautilus_config;

    xml += PairToXML(device2config, "", XML_GNAUTILUS_PAIR_SECOND_CONFIG_META_INFO);
    xml += XML_GNAUTILUS_END_CONFIG_SINGLE_DEVICE;
    xml += GDS_XML_VALUE_NODE_END;
    xml += GDS_XML_METAINFO_END;

    return xml;
}