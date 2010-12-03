#ifndef _dacq_collector_h
#define _dacq_collector_h
#include <string.h>
/*
 * Services
 */
#define DACQ_COLLECTOR   "collector"
#define DACQ_ISOTRAK     "isotrak"
#define DACQ_XPLOTTER    "xplotter"
#define DACQ_JANITOR     "janitor"
#define DACQ_FTP         "dacq-ftp"
/*
 *  Collector commands and response codes
 *
 */
#define COLL_PASSWORD   "pass homunculus122"
#define COLL_SETUP      "setu"
#define COLL_STATUS     "stat"
#define COLL_MEASURE    "meas"
#define COLL_STOP       "stop"
#define COLL_EXIT       "exit"
#define COLL_QUIT       "quit"
#define COLL_SUSPEND    "susp"
#define COLL_RESUME     "resu"
#define COLL_SUSPEND_DISK    "susd"
#define COLL_RESUME_DISK     "resd"
#define COLL_START_AGAIN "sagn"
#define COLL_WAKEUP      "wkup"

#define COLL_STIM_ON    "ston"
#define COLL_STIM_OFF   "stof"

#define COLL_START_PARS "para"
#define COLL_START_STIM "stim"
#define COLL_END_PARS   "."

#define COLL_LAST_ERROR "lerr"

#define COLL_OK               200
#define COLL_SLAVE_NOT_ACTIVE 210
#define COLL_SLAVE_NOT_SET_UP 211
#define COLL_SLAVE_SET_UP     212
#define COLL_MEASURING        213
#define COLL_REPLAYING        214

#define COLL_WARNING(x) ((x) / 100 == 3)

#endif
