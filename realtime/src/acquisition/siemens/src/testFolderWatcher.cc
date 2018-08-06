/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <stdio.h>

#include <FolderWatcher.h>

int main(int argc, char *argv[]) {
	FolderWatcher FW("D:\\watch");
	
	while (true) {
		int num = FW.waitForChanges();
		printf("numChanges = %i\n", num);
		if (num>0) {
			const std::vector<std::string>& vfn = FW.getFilenames();
			printf("VEC[0] = %s\n", vfn[0].c_str());
		}
	}
}
