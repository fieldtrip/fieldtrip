#include <stdlib.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>     

#include <mach/task.h>


void process(struct rusage *p, char *when)
{
		printf("%s\n", when);
		printf(" /* user time used */                   %8ld  %8ld\n",  p->ru_utime.tv_sec,p->ru_utime.tv_usec   );
		printf(" /* system time used */                 %8ld  %8ld\n",  p->ru_stime.tv_sec,p->ru_stime.tv_usec   );
		printf(" /* integral shared memory size */      %8ld\n",  p->ru_ixrss           );
		printf(" /* integral unshared data  */          %8ld\n",  p->ru_idrss           );
		printf(" /* integral unshared stack  */         %8ld\n",  p->ru_isrss           );
		printf(" /* page reclaims */                    %8ld\n",  p->ru_minflt          );
		printf(" /* page faults */                      %8ld\n",  p->ru_majflt          );
		printf(" /* swaps */                            %8ld\n",  p->ru_nswap           );
		printf(" /* block input operations */           %8ld\n",  p->ru_inblock         );
		printf(" /* block output operations */          %8ld\n",  p->ru_oublock         );
		printf(" /* messages sent */                    %8ld\n",  p->ru_msgsnd          );
		printf(" /* messages received */                %8ld\n",  p->ru_msgrcv          );
		printf(" /* signals received */                 %8ld\n",  p->ru_nsignals        );
		printf(" /* voluntary context switches */       %8ld\n",  p->ru_nvcsw           );
		printf(" /* involuntary  */                     %8ld\n",  p->ru_nivcsw          );

}

int getmem (unsigned int *rss, unsigned int *vs) {
		task_t task = MACH_PORT_NULL;
		struct task_basic_info t_info;
		mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

		if (KERN_SUCCESS != task_info(mach_task_self(),
								TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))
		{
				return -1;
		}
		*rss = t_info.resident_size;
		*vs  = t_info.virtual_size;
		return 0;
}


int main()
{
		int ret;
		char *buf;
		int i=0;
		unsigned int rss, vs;
		int who= RUSAGE_SELF;
		struct rusage usage;
		struct rusage *p=&usage;

		while (1) {
				ret=getrusage(who,p);
				/*				process(p, "-------------"); */
				getmem(&rss, &vs);
				printf(" getmem = %lu    %lu\n", rss, vs);
				buf = malloc(1024);
		}

		return 0;
}
