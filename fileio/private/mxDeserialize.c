/*
 * mxDeserialize wrapper
 *
 * Converts a uint8 array to a matlab object, assuming the array has
 * first been created using mxSerialize.
 *
 * Copyright (C) 2005, Brad Phelan         http://xtargets.com
 * Copyright (C) 2007, Robert Oostenveld   http://www.fcdonders.ru.nl
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *
 * $Log: mxDeserialize.c,v $
 * Revision 1.1  2009/01/14 09:43:37  roboos
 * moved source code for mex file from fileio/mex to file/private
 * compiling the source code from within Matlab now ensures that the mex file will be located immediately at the right position
 *
 * Revision 1.2  2007/11/07 12:48:37  roboos
 * added declaration of mx functions, since not included in the matlab header files
 *
 * Revision 1.1  2007/11/07 11:57:55  roboos
 * the original code for this function is released under the LGPL v2
 * Copyright (C) 2005, Brad Phelan, http://xtargets.com
 * renamed the functions to the mxSerialize and mxDeserialize
 * added a check on the number of input and output arguments
 *
 *
 */

#include "mex.h"

/* Only define EXTERN_C if it hasn't been defined already. This allows
 * individual modules to have more control over managing their exports.
 */
#ifndef EXTERN_C
#ifdef __cplusplus
  #define EXTERN_C extern "C"
#else
  #define EXTERN_C extern
#endif
#endif

EXTERN_C mxArray* mxSerialize(const mxArray*);
EXTERN_C mxArray* mxDeserialize(const void*, size_t);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* mxDeserialize is an undocumented Matlab function and should be
    * used assuming the Mathworks may change or remove this function
    * completely from future version of matlab */
  if (nlhs && nrhs)
    plhs[0] = ( mxArray * ) mxDeserialize(mxGetData(prhs[0]), mxGetNumberOfElements(prhs[0]));
}

