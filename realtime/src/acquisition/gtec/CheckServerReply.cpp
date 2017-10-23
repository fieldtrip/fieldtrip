//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#include "CheckServerReply.hpp"

int CheckServerReply(const std::string& xml, bool disp, bool enable_throw) throw (int) 
{
    std::string str_error = ParseXML(xml, GDS_XML_ERROR_NODE);
    std::string status = ParseXML(xml, GDS_XML_STATUS_NODE);

    int error = 0;
    
	if (status != STATUS_SUCCESS || str_error != ERROR_NO_ERROR) 
	{
        error = (str_error.empty()) ? -1 : atoi(str_error.c_str());
        
		if (disp) 
		{
            std::string error_message = ParseXML(EscapeXML(ParseXML(xml, GDS_XML_PAYLOAD_NODE), false), GDS_XML_VALUE_NODE);
            std::cerr << "Server reports an error " << ParseXML(xml, GDS_XML_COMMAND_NODE) << std::endl;
            std::cerr << "Command: " << std::endl;
            std::cerr << "Error message: " << error_message << std::endl;
            std::cerr << "Status code: " << status << std::endl;
            std::cerr << "Error code: " << error << std::endl;
        }

        if (enable_throw)
            throw error;
    }

    return error;
}
