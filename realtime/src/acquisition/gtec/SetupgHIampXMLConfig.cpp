#include "SetupgHIampXMLConfig.hpp"

std::string SetupgHIampXMLConfig(std::string device_name,
        std::string sampling_rate, 
		size_t number_of_channels,
        std::string number_of_scans, 
		std::string trigger_lines_enabled,
        std::string hold_enabled, 
		std::string is_slave,
        std::string counter_enabled, 
		std::string signal_generator_enabled,
        std::string wave_shape, 
		std::string amplitude, 
		std::string frequency,
        std::string offset) 
{
    std::string xml;

    xml = GDS_XML_DECLARATION;
    xml += GDS_XML_DOCTYPE;
    xml += GDS_XML_METAINFO_BEGIN;
    xml += GDS_XML_VALUE_NODE_BEGIN;
    xml += XML_GHIAMP_BEGIN_CONFIG_SINGLE_DEVICE;

    std::string ghiamp_config(XML_GHIAMP_BASE_META_INFO_BEGIN);

    ghiamp_config += XML_GHIAMP_SAMPLE_RATE_BEGIN + sampling_rate + XML_GHIAMP_SAMPLE_RATE_END; //sampling_rate + XML_GHIAMP_SAMPLE_RATE_END;
    ghiamp_config += XML_GHIAMP_NUMBER_OF_SCANS_BEGIN + number_of_scans + XML_GHIAMP_NUMBER_OF_SCANS_END;
    ghiamp_config += XML_GHIAMP_TRIGGER_LINES_ENABLED_BEGIN + trigger_lines_enabled + XML_GHIAMP_TRIGGER_LINES_ENABLED_END;
    ghiamp_config += XML_GHIAMP_HOLD_ENABLED_BEGIN + hold_enabled + XML_GHIAMP_HOLD_ENABLED_END;
    ghiamp_config += XML_GHIAMP_SLAVE_MODE_BEGIN + is_slave + XML_GHIAMP_SLAVE_MODE_END;
    ghiamp_config += XML_GHIAMP_COUNTER_ENABLED_BEGIN + counter_enabled + XML_GHIAMP_COUNTER_ENABLED_END;

    std::string internal_signal_generator(XML_GHIAMP_PARENT_INTERNAL_SIGNAL_GENERATOR_BEGIN);
    internal_signal_generator += XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_ENABLED_BEGIN + signal_generator_enabled + XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_ENABLED_END;
    internal_signal_generator += XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_WAVE_SHAPE_BEGIN + wave_shape + XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_WAVE_SHAPE_END;
    internal_signal_generator += XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_AMPLITUDE_BEGIN + amplitude + XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_AMPLITUDE_END;
    internal_signal_generator += XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_FREQUENCY_BEGIN + frequency + XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_FREQUENCY_END;
    internal_signal_generator += XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_OFFSET_BEGIN + offset + XML_GHIAMP_INTERNAL_SIGNAL_GENERATOR_OFFSET_END;
    internal_signal_generator += XML_GHIAMP_PARENT_INTERANL_SIGNAL_GENERATOR_END;
    ghiamp_config += internal_signal_generator;

    ghiamp_config += XML_GHIAMP_PARENT_CHANNELS_BEGIN;

    std::vector<std::string> channels;
    for (size_t i = 0; i < 256; i++) 
	{
        std::stringstream sstr;
        sstr << i + 1;
        std::string channel_setup;
        channel_setup = XML_GHIAMP_CHANNEL_NUMBER_BEGIN + sstr.str() + XML_GHIAMP_CHANNEL_NUMBER_END;
        sstr.str("");
        sstr << (i < number_of_channels);
        channel_setup += XML_GHIAMP_CHANNEL_ACQUIRE_BEGIN + sstr.str() + XML_GHIAMP_CHANNEL_ACQUIRE_END;
        channel_setup += XML_GHIAMP_CHANNEL_BANDPASS_FILTER_INDEX_BEGIN + std::string("-1") + XML_GHIAMP_CHANNEL_BANDPASS_FILTER_INDEX_END;
        channel_setup += XML_GHIAMP_CHANNEL_NOTCH_FILTER_INDEX_BEGIN + std::string("-1") + XML_GHIAMP_CHANNEL_NOTCH_FILTER_INDEX_END;
        channel_setup += XML_GHIAMP_CHANNEL_REFERENCE_BEGIN + std::string("3") + XML_GHIAMP_CHANNEL_REFERENCE_END;
        channels.push_back(channel_setup);
    }

    ghiamp_config += VectorToXML(channels, std::string(XML_GHIAMP_CHANNEL_VECTOR_ITEM_EXTENSION));

    ghiamp_config += XML_GHIAMP_PARENT_CHANNELS_END;

    ghiamp_config += XML_GHIAMP_BASE_META_INFO_END;

    std::pair<std::string, std::string> device2config;
    device2config.first = device_name;
    device2config.second = ghiamp_config;

    xml += PairToXML(device2config, "", XML_GHIAMP_PAIR_SECOND_CONFIG_META_INFO);
    xml += XML_GHIAMP_END_CONFIG_SINGLE_DEVICE;
    xml += GDS_XML_VALUE_NODE_END;
    xml += GDS_XML_METAINFO_END;

    return xml;
}
