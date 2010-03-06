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

int getmem (unsigned int *rss, unsigned int *vs) {

#if defined (PLATFORM_OSX)

		task_t task = MACH_PORT_NULL;
		struct task_basic_info t_info;
		mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
		if (KERN_SUCCESS != task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count)) {
				return -1;
		}
		*rss = t_info.resident_size;
		*vs  = t_info.virtual_size;
		return 0;

#elif defined(PLATFORM_WIN32)

#elif defined(PLATFORM_LINUX)

#endif

		/* return an error if non of the defines worked out */
		return -1;

}

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
