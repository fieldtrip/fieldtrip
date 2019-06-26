//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#include "EscapeXML.hpp"

std::string EscapeXML(std::string xml, bool escape) 
{
    if (xml.empty()) 
	{
        std::cerr << "ERROR: the provided XML expression is empty" << std::endl;
        return "";
    }

    const static std::string xml_token_quot(XML_TOKEN_QUOT);
    const static std::string xml_token_apos(XML_TOKEN_APOS);
    const static std::string xml_token_lt(XML_TOKEN_LT);
    const static std::string xml_token_gt(XML_TOKEN_GT);
    const static std::string xml_token_amp(XML_TOKEN_AMP);

    const static std::string xml_escape_quot(XML_ESCAPE_QUOT);
    const static std::string xml_escape_apos(XML_ESCAPE_APOS);
    const static std::string xml_escape_lt(XML_ESCAPE_LT);
    const static std::string xml_escape_gt(XML_ESCAPE_GT);
    const static std::string xml_escape_amp(XML_ESCAPE_AMP);

    size_t pos = 0;
    if (escape) 
	{
        pos = xml.find(xml_token_amp);
        while (pos != std::string::npos) 
		{
            xml.replace(pos, xml_token_amp.length(), xml_escape_amp);
            pos = xml.find(xml_token_amp, pos + xml_escape_amp.size());
        }

        pos = xml.find(xml_token_quot);
        while (pos != std::string::npos) 
		{
            xml.replace(pos, xml_token_quot.length(), xml_escape_quot);
            pos = xml.find(xml_token_quot, pos + xml_escape_quot.size());
        }

        pos = xml.find(xml_token_apos);
        while (pos != std::string::npos) 
		{
            xml.replace(pos, xml_token_apos.length(), xml_escape_apos);
            pos = xml.find(xml_token_apos, pos + xml_escape_apos.size());
        }

        pos = xml.find(xml_token_lt);
        while (pos != std::string::npos) 
		{
            xml.replace(pos, xml_token_lt.length(), xml_escape_lt);
            pos = xml.find(xml_token_lt, pos + xml_escape_lt.size());
        }

        pos = xml.find(xml_token_gt);
        while (pos != std::string::npos) 
		{
            xml.replace(pos, xml_token_gt.length(), xml_escape_gt);
            pos = xml.find(xml_token_gt, pos + xml_escape_gt.size());
        }

        return xml;
    }

    // unescape XML
    pos = xml.find(xml_escape_amp);
    while (pos != std::string::npos) 
	{
        xml.replace(pos, xml_escape_amp.length(), xml_token_amp);
        pos = xml.find(xml_escape_amp, pos + xml_token_amp.size());
    }

    pos = xml.find(xml_escape_quot);
    while (pos != std::string::npos) 
	{
        xml.replace(pos, xml_escape_quot.length(), xml_token_quot);
        pos = xml.find(xml_escape_quot, pos + xml_token_quot.size());
    }

    pos = xml.find(xml_escape_apos);
    while (pos != std::string::npos) 
	{
        xml.replace(pos, xml_escape_apos.length(), xml_token_apos);
        pos = xml.find(xml_escape_apos, pos + xml_token_apos.size());
    }

    pos = xml.find(xml_escape_lt);
    while (pos != std::string::npos) 
	{
        xml.replace(pos, xml_escape_lt.length(), xml_token_lt);
        pos = xml.find(xml_escape_lt, pos + xml_token_lt.size());
    }

    pos = xml.find(xml_escape_gt);
    while (pos != std::string::npos) 
	{
        xml.replace(pos, xml_escape_gt.length(), xml_token_gt);
        pos = xml.find(xml_escape_gt, pos + xml_token_gt.size());
    }

    return xml;
}
