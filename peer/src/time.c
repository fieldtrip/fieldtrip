/*
 * Copyright (C) 2010, Robert Oostenveld
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include "mex.h"


void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {

		/* return the time since the Epoch (00:00:00 UTC, January 1, 1970), measured in seconds */
		plhs[0] = mxCreateDoubleScalar(time(NULL));

		return;
} /* main */

