#include <syslog.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

		printf("LOG_EMERG      %d\n", LOG_EMERG    );
		printf("LOG_ALERT      %d\n", LOG_ALERT    );
		printf("LOG_CRIT       %d\n", LOG_CRIT     );
		printf("LOG_ERR        %d\n", LOG_ERR      );
		printf("LOG_WARNING    %d\n", LOG_WARNING  );
		printf("LOG_NOTICE     %d\n", LOG_NOTICE   );
		printf("LOG_INFO       %d\n", LOG_INFO     );
		printf("LOG_DEBUG      %d\n", LOG_DEBUG    );

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
