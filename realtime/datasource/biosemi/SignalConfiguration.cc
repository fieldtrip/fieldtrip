#include <SignalConfiguration.h>
#include <stdio.h>
#include <string.h>

int SignalConfiguration::parseFile(const char *filename) {
	FILE *fp;
	char line[2048];
	int lineCount = 0;
	int errorCount = 0;
	int addSave = 1, addStream = 1;

	chanSelSave.clear();
	chanSelStream.clear();
	maxChanSave = maxChanStream = 0;

	fp = fopen(filename, "r");
	if (fp==NULL) return -1;
	
	while (!feof(fp)) {
		if (fgets(line, sizeof(line), fp) == NULL) break;
		lineCount++;
	
		char *lp = line;
		while (isspace(*lp)) lp++;
		// ignore comments starting by ; or # as well as empty lines
		if (*lp == ';' || *lp == '#' || *lp==0) continue;
		if (strchr(lp, '=')) {
			int n,chn; 
			char lab[100];
			n = sscanf(lp, "%i=%99s", &chn, lab);
			if (n==2 && chn>0) {
				if (addSave) {
					if (chn>maxChanSave) maxChanSave = chn;
					chanSelSave.add(chn-1, lab);
				}
				if (addStream) {
					if (chn>maxChanStream) maxChanStream = chn;
					chanSelStream.add(chn-1, lab);
				}
				continue;
			}
		}
		if (!strncasecmp(lp, "downsample", 10)) {
			int ds = atoi(lp+10);
			if (ds>=1) {
				downSample=ds;
				continue;
			}
		}
		if (!strncasecmp(lp, "[select]", 8)) {
			addStream = addSave = 1;
			continue;
		}
		if (!strncasecmp(lp, "[stream]", 8)) {
			addStream = 1;
			addSave = 0;
			continue;
		}
		if (!strncasecmp(lp, "[save]", 6)) {
			addSave = 1;
			addStream = 0;
			continue;
		}
		fprintf(stderr, "Could not parse line %i:\n%s\n", lineCount, line);
		errorCount++;
	}
	fclose(fp);
	return errorCount;
}
