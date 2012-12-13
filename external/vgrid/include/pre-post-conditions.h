#ifndef NMWR_GB_PRE_POST_CONDITION_H
#define NMWR_GB_PRE_POST_CONDITION_H

/*
  © Copyright 2003, C&C Research Laboratories, NEC Europe Ltd.
  

    This file is part of SimBio-Vgrid.

    SimBio-Vgrid is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SimBio-Vgrid is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SimBio-Vgrid.  If not, see <http://www.gnu.org/licenses/>.


*/


#include <iostream>
#include <stdlib.h>

  
//----------------------------------------------------------------
/*! \file 
 

   \brief Some useful macros for checking pre/post conditions.

   $Id$


  Two different sets of macros are defined:
 - REQUIRE/ENSURE (condition, message, severity)
   are defined only if  \c NMWR_DEBUG is defined
 - REQUIRE_ALWAY/ENSURE_ALWAYS (condition, message, severity)
   are always defined

 Example:
 \code
 int* f(int* p) {
  REQUIRE((p != 0), "Pointer null!", 1);
  // work on p ...
  ENSURE((p != 0), "Pointer null!", 1);
 }
 \endcode
*/ 
//----------------------------------------------------------------


#define _ERRORLOG std::cerr

#define _PRECONDITION_ERROR  _ERRORLOG << "\nERROR in FILE "  << __FILE__ << ", LINE " << __LINE__\
   << "\n(compiled on " << __DATE__ << " at " __TIME__  << " )\n"\
   << "precondition violated:\n"

#define _POSTCONDITION_ERROR  _ERRORLOG << "\nERROR in FILE "  << __FILE__ << ", LINE " << __LINE__\
   << "\n(compiled on " << __DATE__ << " at " __TIME__  << " )\n"\
                                       << ": postcondition violated:\n"

#define REQUIRE_ALWAYS(condition, error_msg, severity)\
 if(! (condition))  { _PRECONDITION_ERROR << #condition << ' ' << error_msg << std::endl; abort();}

#define ENSURE_ALWAYS(condition, error_msg, severity)\
 if(! (condition)) { _POSTCONDITION_ERROR << #condition << ' ' << error_msg << std::endl; abort();}

#ifndef NDEBUG

#define REQUIRE(condition, error_msg, severity)\
 if(! (condition))  { _PRECONDITION_ERROR << #condition << ' ' << error_msg << std::endl; abort();}

#define ENSURE(condition, error_msg, severity)\
 if(! (condition)) { _POSTCONDITION_ERROR << #condition << ' '  << error_msg << std::endl; abort();}

#else
#define REQUIRE(condition, error_msg, severity) 
#define ENSURE(condition, error_msg, severity)
#endif



#endif
