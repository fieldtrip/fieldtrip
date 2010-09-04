#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <pwd.h>

int main()
{
		struct passwd *pwd;

		pwd = getpwuid(geteuid());
		printf("user = %s\n", pwd->pw_name);

		return 0;
}
