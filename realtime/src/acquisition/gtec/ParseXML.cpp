//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#include "ParseXML.hpp"

std::string ParseXML(std::string xml, std::string xml_node, size_t* pos_finally, size_t pos_begin = 0) 
{
    if (xml.find(GDS_XML_DOCTYPE) == std::string::npos)
	{
        std::cerr << "ERROR: did not find " << GDS_XML_DOCTYPE<< std::endl;
        return "";
    }

    std::string begin_token(XML_TOKEN_LT + xml_node + XML_TOKEN_GT);
    std::string end_token(XML_TOKEN_LTS + xml_node + XML_TOKEN_GT);

    pos_begin = xml.find(begin_token, pos_begin);
    if (pos_begin == std::string::npos) 
	{
        std::cerr << "ERROR: could not find the token " << begin_token << std::endl;
        return "";
    }

    size_t pos_end = xml.find(end_token, pos_begin + begin_token.size());
    if (pos_end == std::string::npos) 
	{
        std::cerr << "ERROR: could not find the token " << end_token << std::endl;
        return "";
    }

    std::string ret(xml.begin() + pos_begin + begin_token.size(), xml.begin() + pos_end);
    *pos_finally = pos_end + end_token.size();

    return ret;
}

std::string ParseXML(std::string xml, std::string xml_node) 
{
    size_t dummy = 0;
    return ParseXML(xml, xml_node, &dummy, dummy);
}

void ParseXML(std::string xml, size_t* scan_count, size_t* channel_count, size_t* buffer_size_in_samples) 
{
    std::string data_info = EscapeXML(ParseXML(xml, GDS_XML_PAYLOAD_NODE), false);
    std::vector<std::string> vct_channles_per_device = XMLToVector(ParseXML(data_info, GDS_XML_DATA_INFO_CHANNELS_PER_DEVICE_NODE));
    *channel_count = atoi(vct_channles_per_device[0].c_str());
    *scan_count = atoi(ParseXML(data_info, GDS_XML_DATA_INFO_NUMBER_OF_SCANS_NODE).c_str());
    *buffer_size_in_samples = atoi(ParseXML(data_info, GDS_XML_DATA_INFO_BUFFER_SIZE_IN_SAMPLES_NODE).c_str());
}
