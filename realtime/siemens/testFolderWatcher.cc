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
