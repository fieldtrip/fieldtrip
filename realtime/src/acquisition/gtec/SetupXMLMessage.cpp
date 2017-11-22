//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#include "SetupXMLMessage.hpp"

std::string SetupXMLElement(std::string xml_node, std::string value, bool new_line)
{
    std::string token_close = new_line ? ">\n" : ">";
    return "<" + xml_node + token_close + value + "</" + xml_node + ">\n";
}

std::string SetupXMLElement(std::string xml_node, size_t value, bool new_line)
{
    std::stringstream to_str;
    to_str << value;
    std::string token_close = new_line ? ">\n" : ">";
    return "<" + xml_node + token_close + to_str.str() + "</" + xml_node + ">\n";
}

std::string SetupXMLMessage(std::string session_id, std::string command, std::string payload)
{
    std::string message;

    message = GDS_XML_DECLARATION;
    message += GDS_XML_DOCTYPE;
    message += GDS_XML_METAINFO_BEGIN;
    message += GDS_XML_VALUE_NODE_BEGIN;
    message += SetupXMLElement(GDS_XML_SESSION_ID_NODE, session_id);
    message += SetupXMLElement(GDS_XML_COMMAND_NODE, command);

    if (payload.empty())
        message += SetupXMLElement(GDS_XML_PAYLOAD_NODE, "");
    else
        message += SetupXMLElement(GDS_XML_PAYLOAD_NODE, payload);

    message += SetupXMLElement(GDS_XML_STATUS_NODE, STATUS_COMMAND);
    message += SetupXMLElement(GDS_XML_ERROR_NODE, ERROR_NO_ERROR);
    message += GDS_XML_VALUE_NODE_END;
    message += GDS_XML_METAINFO_END;

    return message;
}

std::string SetupXMLMessage(std::string xml_payload)
{
    std::string xml;

    xml = GDS_XML_DECLARATION;
    xml += GDS_XML_DOCTYPE;
    xml += GDS_XML_METAINFO_BEGIN;
    xml += GDS_XML_VALUE_NODE_BEGIN;
    xml += xml_payload;
    xml += GDS_XML_VALUE_NODE_END;
    xml += GDS_XML_METAINFO_END;

    return xml;
}

std::string SetupXMLMessage(size_t scan_count, size_t channel_count, size_t buffer_size_in_samples)
{
    std::stringstream to_str;
    to_str << channel_count;
    std::string channels_per_device = (channel_count == 0) ? VectorToXML(std::vector<std::string>()) : VectorToXML(std::vector<std::string>(1, to_str.str()));

    std::string xml_payload;
    xml_payload = SetupXMLElement("number_of_scans_", scan_count);
    xml_payload += SetupXMLElement("channels_per_device_", channels_per_device, true);
    xml_payload += SetupXMLElement("buffer_size_in_samples_", buffer_size_in_samples);

    return SetupXMLMessage(xml_payload);
}
