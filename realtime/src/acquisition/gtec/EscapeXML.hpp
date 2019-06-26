//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef ESCAPEXML_HPP_INCLUDED
#define ESCAPEXML_HPP_INCLUDED

#include <string>
#include <iostream>
#include "XMLDefinitions.hpp"

// Replaces XML specific charecters with their escape code.
// If escape is set to false, then the xml is interpreted as escaped XML and
// the escape code is replaced by the XML specific character.
//-------------------------------------------------------------------------------------
std::string EscapeXML(std::string xml, bool escape = true);

#endif // ESCAPEXML_HPP_INCLUDED
