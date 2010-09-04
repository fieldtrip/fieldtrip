#include <syslog.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

		openlog("test", LOG_PERROR | LOG_PID, LOG_USER);
		syslog(LOG_EMERG        , "LOG_EMERG    ");
		syslog(LOG_ALERT        , "LOG_ALERT    ");
		syslog(LOG_CRIT         , "LOG_CRIT     ");
		syslog(LOG_ERR          , "LOG_ERR      ");
		syslog(LOG_WARNING      , "LOG_WARNING  ");
		syslog(LOG_NOTICE       , "LOG_NOTICE   ");
		syslog(LOG_INFO         , "LOG_INFO     ");
		syslog(LOG_DEBUG        , "LOG_DEBUG    ");
		closelog();


		return 0;
}
