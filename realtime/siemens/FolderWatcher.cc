#include <FolderWatcher.h>

FolderWatcher::FolderWatcher(const char *directory) : vecFilenames(0) {
	dirHandle = INVALID_HANDLE_VALUE;
	isListening = false;
	
	eventHandle = CreateEvent(NULL, FALSE, FALSE, NULL);
	if (eventHandle == INVALID_HANDLE_VALUE) return;
	
	dirHandle = CreateFile(directory, 
					FILE_LIST_DIRECTORY, 
					FILE_SHARE_WRITE|FILE_SHARE_READ|FILE_SHARE_DELETE, 
					NULL, 
					OPEN_EXISTING, 
					FILE_FLAG_BACKUP_SEMANTICS|FILE_FLAG_OVERLAPPED, 
					NULL);
}

int FolderWatcher::waitForChanges() {
	DWORD bufferLength=0;
	
	vecFilenames.clear();
	
	if (dirHandle==INVALID_HANDLE_VALUE) return -1;
		
	BOOL ok = ReadDirectoryChangesW(dirHandle, fileInfoBuffer, sizeof(fileInfoBuffer), TRUE, FILE_NOTIFY_CHANGE_SIZE, &bufferLength, NULL, NULL);
	
	if (!ok) return 0;	
	if (bufferLength==0) return 0;
	
	return processChanges();
}


bool FolderWatcher::startListenForChanges(HANDLE extHandle) {
	vecFilenames.clear();
	
	if (dirHandle==INVALID_HANDLE_VALUE) return 0;
	
	if (extHandle==INVALID_HANDLE_VALUE) {
		usedHandle = overlap.hEvent = eventHandle;
	} else {
		usedHandle = overlap.hEvent = extHandle;
	}
	
	if (ReadDirectoryChangesW(dirHandle, fileInfoBuffer, sizeof(fileInfoBuffer), TRUE, FILE_NOTIFY_CHANGE_LAST_WRITE, NULL, &overlap, NULL)) {
		isListening = true;
		return true;
	}
	return false;
}

bool FolderWatcher::stopListenForChanges() {
	if (isListening) {
		if (CancelIo(dirHandle)) {
			isListening = false;
			return true;
		}
		return false;
	}
	return true;
}

int FolderWatcher::checkHasChanged(unsigned int millis) {
	if (WaitForSingleObject(usedHandle, millis) != WAIT_OBJECT_0) return 0;
	
	return processChanges();
}
	
int FolderWatcher::processChanges() {
	FILE_NOTIFY_INFORMATION *info;

	char *fileInfoBufferPtr = fileInfoBuffer;
	
	memset(fileInfoBuffer, 0, sizeof(DWORD));
	isListening = false;	
	
	do {
		char filename[MAX_PATH];
		info = (FILE_NOTIFY_INFORMATION *) fileInfoBufferPtr;
		
		fileInfoBufferPtr += info->NextEntryOffset;
		/*
		switch(info->Action) {
			case FILE_ACTION_ADDED:
				printf("File added\n");
				break;
			case FILE_ACTION_MODIFIED:
				printf("File modified\n");
				break;
			default:
				printf("Other action\n");
		}
		*/
		#ifdef UNICODE
        lstrcpynW(filename, info->FileName, min(MAX_PATH, info->FileNameLength / sizeof(WCHAR) + 1));
		#else
        {
            int count = WideCharToMultiByte(CP_ACP, 0, info->FileName, info->FileNameLength / sizeof(WCHAR), filename, MAX_PATH - 1, NULL, NULL);
            filename[count] = 0;
        }
		#endif
		vecFilenames.push_back(filename);
		printf("%i: %s\n", vecFilenames.size(), filename);
	} while (info->NextEntryOffset > 0);
	
	return vecFilenames.size();
}

FolderWatcher::~FolderWatcher() {	
	if (isListening) {
		CancelIo(dirHandle);
	}
	CloseHandle(dirHandle);
	CloseHandle(eventHandle);
}
