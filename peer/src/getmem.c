/*
 * Copyright (C) 2010, Robert Oostenveld
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/
 *
 */

#include <stdlib.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>     
#include <mach/task.h>

#include "platform.h"
#include "mex.h"
#include "matrix.h"

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
		unsigned int rss, vs;

		if (getmem(&rss, &vs)==0) {
				/* mexPrintf(" getmem = %lu    %lu\n", rss, vs); */
				plhs[0] = mxCreateDoubleScalar(rss);
				plhs[1] = mxCreateDoubleScalar(vs);
		}
		else {
				plhs[0] = mxCreateDoubleScalar(mxGetNaN());
				plhs[1] = mxCreateDoubleScalar(mxGetNaN());
		}

		return;
}
