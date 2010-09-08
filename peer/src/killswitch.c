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

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

int check_killswitch(void) {
		int found = 0;
		peerlist_t *peer;

		pthread_mutex_lock(&mutexkillswitch);
		if (killswitch.enabled==0) {
				pthread_mutex_unlock(&mutexkillswitch);
				return 0;
		}
		if (killswitch.masterid==0) {
				pthread_mutex_unlock(&mutexkillswitch);
				return 0;
		}

		pthread_mutex_lock(&mutexpeerlist);
		peer = peerlist;
		while(peer && !found) {
				found = 1;
				found = found && (peer->host->id == killswitch.masterid);
				found = found && (peer->host->status == STATUS_MASTER);
				peer = peer->next;
		}

		pthread_mutex_unlock(&mutexpeerlist);
		pthread_mutex_unlock(&mutexkillswitch);

		if (!found) {
				DEBUG(LOG_NOTICE, "the kill switch was triggered");
				exit(1);
		}

		return 0;
}

