//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#include "XMLSTLConversions.hpp"

std::vector<std::string> XMLToVector(std::string xml)
{
    size_t dummy = 0;
    return XMLToVector(xml, &dummy, 0);
}

std::vector<std::string> XMLToVector(std::string xml, size_t* pos_finally, size_t pos_init)
{
    std::vector<std::string> ret;
    *pos_finally = 0;

    size_t pos_begin = xml.find(token_count_begin, pos_init);
	if (pos_begin == std::string::npos)
        pos_begin = xml.find(token_count_begin, pos_init);

    if (pos_begin == std::string::npos)
	{
        std::cerr << "ERROR: xml expression does not hold a list" << std::endl;
        return ret;
    }

    size_t pos_end = xml.find(token_count_end, pos_begin + size_token_count_begin);
    size_t count = atoi(std::string(xml.begin() + pos_begin + size_token_count_begin, xml.begin() + pos_end).c_str());

    pos_begin = xml.find(xml_stl_list_item_version, pos_end + token_count_end.size());
    if (pos_begin == std::string::npos)
	{
        std::cerr << "ERROR: xml expression does not hold a list" << std::endl;
        return ret;
    }

    pos_begin = xml.find(token_item_begin, pos_begin + xml_stl_list_item_version.size());
    pos_end = xml.find(token_item_end, pos_begin + size_token_begin);
    size_t index = 0;
    while (pos_begin != std::string::npos)
	{
        index++;
        ret.push_back(std::string(xml.begin() + pos_begin + size_token_begin, xml.begin() + pos_end));

        if (index == count)
            break;
        
		pos_begin = xml.find(token_item_begin, pos_end + size_token_end);
        pos_end = xml.find(token_item_end, pos_begin + size_token_begin);
    }

    if (index != count)
        std::cerr << "ERROR: Expected " << count << " number of items, but found " << index << " number of items." << std::endl;

    *pos_finally = pos_end + size_token_end;

    return ret;
}

std::vector<struct Device> XMLToDevice(std::string xml)
{
    std::vector<struct Device> device_list;
    std::string inner_count;
    std::string status;
    std::string name;
    std::string type;
    size_t pos_begin = 0;
    size_t pos_init = 0;
    size_t pos_end = 0;

    std::string count = ParseXML(xml, XML_STL_LIST_COUNT, &pos_end, pos_begin);
    for (int i = 0; i < atoi(count.c_str()); i++)
	{
        status = ParseXML(xml, GDS_XML_DEVICE_INFO_STATUS, &pos_end, pos_init);
        inner_count = ParseXML(xml, XML_STL_LIST_COUNT, &pos_end, pos_end);

        for (int j = 0; j < atoi(inner_count.c_str()); j++)
		{
            name = ParseXML(xml, GDS_XML_DEVICE_INFO_NAME, &pos_end, pos_end);
            type = ParseXML(xml, GDS_XML_DEVICE_INFO_TYPE, &pos_end, pos_end);
            device_list.push_back(Device(name, (GDEVICE_TYPE) atoi(type.c_str()), atoi(status.c_str())));
        }

        pos_init = pos_end;
    }

    return device_list;
}

std::string VectorToXML(std::vector<std::string> v, bool enable_item_version, std::string item_version)
{
    size_t count = v.size();
    std::stringstream sstr_count;
    sstr_count << count;
    std::string xml(token_count_begin + sstr_count.str() + token_count_end);

    if (enable_item_version)
	{
        static const std::string xml_item_version(XML_TOKEN_LT + xml_stl_list_item_version + XML_TOKEN_GT + item_version + XML_TOKEN_LTS + xml_stl_list_item_version + XML_TOKEN_GT + "\n");
        xml += xml_item_version;
    }

    if (v.empty())
        return xml;

    std::vector<std::string>::iterator it = v.begin();
    std::vector<std::string>::iterator end = v.end();

    for (; it != end; ++it)
        xml += token_item_begin + *it + token_item_end;

    return xml;
}

std::string VectorToXML(std::vector<std::string> v, std::string item_extension)
{
    size_t count = v.size();
    std::stringstream sstr_count;
    sstr_count << count;
    std::string xml(token_count_begin + sstr_count.str() + token_count_end);

    if (v.empty())
        return xml;

    std::vector<std::string>::iterator it = v.begin();
    std::vector<std::string>::iterator end = v.end();

    if (!item_extension.empty())
        xml += item_extension + *it + token_item_end;

    ++it;
    for (; it != end; ++it)
        xml += token_item_begin + *it + token_item_end;

    return xml;
}

std::pair<std::string, std::string> XMLToPair(std::string xml)
{
    size_t pos_begin = xml.find(token_first_begin);
    std::pair<std::string, std::string> ret;

    if (pos_begin == std::string::npos)
	{
        std::cerr << "ERROR: xml expression does not contain \"first\" item of std::pair" << std::endl;
        return ret;
    }

    size_t pos_end = xml.find(token_first_end, pos_begin + size_token_first_begin);
    ret.first.append(xml.begin() + pos_begin + size_token_first_begin, xml.begin() + pos_end);
    pos_begin = xml.find(token_second_begin, pos_end + token_first_end.size());

    std::cout << "begin second = " << pos_begin + size_token_second_begin << std::endl;

    if (pos_begin == std::string::npos)
	{
        std::cerr << "ERROR: xml expression does not contain \"second\" item of std::pair" << std::endl;
        return ret;
    }

    pos_end = xml.find(token_second_end, pos_begin + size_token_second_begin);

    std::cout << "end second = " << pos_end << std::endl;

    ret.second.append(xml.begin() + pos_begin + size_token_second_begin, xml.begin() + pos_end);

    return ret;
}

std::string PairToXML(std::pair<std::string, std::string> p, std::string xml_meta_info_first, std::string xml_meta_info_second) {
    std::string xml;

    if (xml_meta_info_first.empty())
        xml += token_first_begin + p.first + token_first_end;
    else 
	{
        std::string token_first_meta_begin(XML_TOKEN_LT + xml_stl_pair_first + xml_meta_info_first + XML_TOKEN_GT + "\n");
        xml += token_first_meta_begin + p.first + token_first_end;
    }

    if (xml_meta_info_second.empty())
        xml += token_second_begin + p.second + token_second_end;
    else 
	{
        std::string token_second_meta_begin(XML_TOKEN_LT + xml_stl_pair_second + xml_meta_info_second + XML_TOKEN_GT + "\n");
        xml += token_second_meta_begin + p.second + token_second_end;
    }

    return xml;
}

