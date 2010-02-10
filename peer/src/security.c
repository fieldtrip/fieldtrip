#include <stdio.h>
#include <stdlib.h>

#include "peer.h"
#include "extern.h"

int security_check(hostdef_t *host) {
		int ismember = 1;

		ismember = ismember & ismember_userlist(host->user);
		ismember = ismember & ismember_grouplist(host->group);
		ismember = ismember & ismember_hostlist(host->name);

		return ismember;
}

int ismember_userlist(char *str) {
  int ismember = 1;
  userlist_t *user = NULL;
  pthread_mutex_lock(&mutexuserlist);
  if (userlist) {
    ismember = 0;
    user = userlist;
    while (user) {
      if (strncmp(user->name, str, STRLEN)==0)
        ismember = 1;
      user = user->next;
    }
  }
  pthread_mutex_unlock(&mutexuserlist);
		return ismember;
}

int ismember_grouplist(char *str) {
		int ismember = 1; 
		grouplist_t *group = NULL;
		pthread_mutex_lock(&mutexgrouplist);
		if (grouplist) {
				ismember = 0;
				group = grouplist;
				while (group) {
						if (strncmp(group->name, str, STRLEN)==0)
								ismember = 1;
						group = group->next;
				}
		}
		pthread_mutex_unlock(&mutexgrouplist);
		return ismember;
}

int ismember_hostlist(char *str) {
		int ismember = 1; 
		hostlist_t *host = NULL;
		pthread_mutex_lock(&mutexhostlist);
		if (hostlist) {
				ismember = 0;
				host = hostlist;
				while (host) {
						if (strncmp(host->name, str, STRLEN)==0)
								ismember = 1;
						host = host->next;
				}
		}
		pthread_mutex_unlock(&mutexhostlist);
		return ismember;
}
