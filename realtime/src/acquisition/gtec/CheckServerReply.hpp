//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef CHECKSERVERREPLY_HPP_INCLUDED
#define CHECKSERVERREPLY_HPP_INCLUDED

#include "ParseXML.hpp"
#include "Constants.hpp"

// Validate the server's reply message.
// If status is not set to succes or an error value is set,
// then retrieve the error code from the data and return
// or throw the error value. Optional further information
// can be displayed.
// If not disabled, then this method throws the error number on error, else
// this method returns the error number.
//-------------------------------------------------------------------------------------
int CheckServerReply(const std::string& xml, bool disp = true, bool enable_throw = true) throw (int);

#endif // CHECKSERVERREPLY_HPP_INCLUDED
