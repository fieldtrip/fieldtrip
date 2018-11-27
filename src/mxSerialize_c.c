/* 
 * mxSerialize wrapper
 *
 * Converts any matlab object into a uint8 array suitable for passing 
 * down a comms channel to be reconstructed at the other end.
 *
 * Copyright (C) 2005, Brad Phelan         http://xtargets.com
 * Copyright (C) 2007, Robert Oostenveld   http://robertoostenveld.nl
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
 * $Id$
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
   /* mxSerialize is an undocumented Matlab function and should be
    * used assuming the Mathworks may change or remove this function
    * completely from future version of matlab */
  if (nlhs && nrhs)
    plhs[0] = (mxArray *) mxSerialize(prhs[0]);
}

