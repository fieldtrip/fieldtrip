//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#include "SetupgUSBampXMLConfig.hpp"

std::string SetupgUSBampXMLConfig(std::string device_name,
        std::string sample_rate, 
		std::string number_of_scans, 
		std::string slave,
        std::string short_cut, 
		std::string device_mode, 
		std::string scan_dio,
        std::string analog_out_shape, 
		std::string analog_out_frequency,
        std::string analog_out_amplitude, 
		std::string analog_out_offset)
{
    std::string xml;

    xml = GDS_XML_DECLARATION;
    xml += GDS_XML_DOCTYPE;
    xml += GDS_XML_METAINFO_BEGIN;
    xml += GDS_XML_VALUE_NODE_BEGIN;
    xml += XML_GUSBAMP_BEGIN_CONFIG_SINGLE_DEVICE;

    std::string gusbamp_config(XML_GUSBAMP_BASE_META_INFO_BEGIN);

    std::string analog_out_config(XML_GUSBAMP_PARENT_ANALOG_OUT_BEGIN);
    analog_out_config += XML_GUSBAMP_ANALOG_OUT_SHAPE_BEGIN + analog_out_shape + XML_GUSBAMP_ANALOG_OUT_SHAPE_END;
    analog_out_config += XML_GUSBAMP_ANALOG_OUT_FREQUENCY_BEGIN + analog_out_frequency + XML_GUSBAMP_ANALOG_OUT_FREQUENCY_END;
    analog_out_config += XML_GUSBAMP_ANALOG_OUT_AMPLITUDE_BEGIN + analog_out_amplitude + XML_GUSBAMP_ANALOG_OUT_AMPLITUDE_END;
    analog_out_config += XML_GUSBAMP_ANALOG_OUT_OFFSET_BEGIN + analog_out_offset + XML_GUSBAMP_ANALOG_OUT_OFFSET_END;
    analog_out_config += XML_GUSBAMP_PARENT_ANALOG_OUT_END;
    gusbamp_config += analog_out_config;

    gusbamp_config += XML_GUSBAMP_SAMPLE_RATE_BEGIN + sample_rate + XML_GUSBAMP_SAMPLE_RATE_END;
    gusbamp_config += XML_GUSBAMP_NUMBER_OF_SCANS_BEGIN + number_of_scans + XML_GUSBAMP_NUMBER_OF_SCANS_END;
    gusbamp_config += XML_GUSBAMP_SLAVE_MODE_BEGIN + slave + XML_GUSBAMP_SLAVE_MODE_END;
    gusbamp_config += XML_GUSBAMP_SHORT_CUT_BEGIN + short_cut + XML_GUSBAMP_SHORT_CUT_END;

    std::string common_ground_config(XML_GUSBAMP_PARENT_COMMON_GROUND_BEGIN);
    common_ground_config += XML_GUSBAMP_VECTOR_DEFINITION_BEGIN;
    std::vector<std::string> common_ground_vector(4, "0");
    common_ground_config += VectorToXML(common_ground_vector, false);
    common_ground_config += XML_GUSBAMP_VECTOR_DEFINITION_END;
    common_ground_config += XML_GUSBAMP_PARENT_COMMON_GROUND_END;
    gusbamp_config += common_ground_config;

    std::string common_reference_config(XML_GUSBAMP_PARENT_COMMON_REFERENCE_BEGIN);
    common_reference_config += XML_GUSBAMP_VECTOR_DEFINITION_BEGIN;
    std::vector<std::string> common_reference_vector(4, "0");
    common_reference_config += VectorToXML(common_reference_vector, false);
    common_reference_config += XML_GUSBAMP_VECTOR_DEFINITION_END;
    common_reference_config += XML_GUSBAMP_PARENT_COMMON_REFERENCE_END;
    gusbamp_config += common_reference_config;

    gusbamp_config += XML_GUSBAMP_DEVICE_MODE_BEGIN + device_mode + XML_GUSBAMP_DEVICE_MODE_END;
    gusbamp_config += XML_GUSBAMP_SCAN_DIO_BEGIN + scan_dio + XML_GUSBAMP_SCAN_DIO_END;

    std::vector<std::string> auto_filter_vector(16, "-1");
    std::string auto_filter = VectorToXML(auto_filter_vector);

    std::string bandpass_config;
    bandpass_config = XML_GUSBAMP_PARENT_BANDPASS_BEGIN;
    bandpass_config += auto_filter;
    bandpass_config += XML_GUSBAMP_PARENT_BANDPASS_END;
    gusbamp_config += bandpass_config;

    std::string notch_config;
    notch_config = XML_GUSBAMP_PARENT_NOTCH_BEGIN;
    notch_config += auto_filter;
    notch_config += XML_GUSBAMP_PARENT_NOTCH_END;
    gusbamp_config += notch_config;

    std::string bipolar_derivation_config(XML_GUSBAMP_PARENT_BIPOLAR_BEGIN);
    std::vector<std::string> bipolar_derivation_vector(16, "0");
    bipolar_derivation_config += VectorToXML(bipolar_derivation_vector);
    bipolar_derivation_config += XML_GUSBAMP_PARENT_BIPOLAR_END;
    gusbamp_config += bipolar_derivation_config;

    std::string analog_in_config(XML_GUSBAMP_PARENT_ANALOG_IN_BEGIN);
    analog_in_config += XML_GUSBAMP_VECTOR_DEFINITION_BEGIN;
    std::vector<std::string> analog_in_vector;
    int max_channels = (scan_dio == "1") ? 17 : 16;
    for (int i = 1; i <= max_channels; i++)
	{
        std::stringstream sstr;
        sstr << i;
        analog_in_vector.push_back(sstr.str());
    }
    analog_in_config += VectorToXML(analog_in_vector);
    analog_in_config += XML_GUSBAMP_VECTOR_DEFINITION_END;
    analog_in_config += XML_GUSBAMP_PARENT_ANALOG_IN_END;
    gusbamp_config += analog_in_config;

    gusbamp_config += XML_GUSBAMP_BASE_META_INFO_END;

    std::pair<std::string, std::string> device2config;
    device2config.first = device_name;
    device2config.second = gusbamp_config;

    xml += PairToXML(device2config, "", XML_GUSBAMP_PAIR_SECOND_CONFIG_META_INFO);
    xml += XML_GUSBAMP_END_CONFIG_SINGLE_DEVICE;
    xml += GDS_XML_VALUE_NODE_END;
    xml += GDS_XML_METAINFO_END;

    return xml;
}
