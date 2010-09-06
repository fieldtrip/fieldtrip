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

#include <stdio.h>
#include <stdlib.h>

#include "peer.h"
#include "extern.h"

/* returns 1 if the host is in the userlist, grouplist and hostlist, 0 if not */
int security_check(hostdef_t *host) {
		int ismember = 1;
		ismember = ismember && ismember_userlist(host->user);
		ismember = ismember && ismember_grouplist(host->group);
		ismember = ismember && ismember_hostlist(host->name);
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
		if (hostlist) {
				ismember = 0;
				host = hostlist;
				while (host) {
						if (strncmp(host->name, str, STRLEN)==0)
								ismember = 1;
						host = host->next;
				}
		}
		return ismember;
}

