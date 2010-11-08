#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>

#include "engine.h"
#include "matrix.h"

int main(int argc, char *argv[]) {
		Engine *en;
		en = engOpen("matlab -nojvm");
		engEvalString(en, "pause(10)");
		engClose(en);

		return 0;
}
