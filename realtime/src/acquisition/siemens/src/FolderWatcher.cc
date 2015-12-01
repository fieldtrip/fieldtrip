/*
 * Copyright (C) 2010, Stefan Klanke
 * 	Modified by Tim van Mourik 2015
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <FolderWatcher.h>

FolderWatcher::FolderWatcher(const char *directory) : vecFilenames(0) {
	dirHandle = INVALID_HANDLE_VALUE;
	isListening = false;
    completionPort = 0;

	activeBuffer = 0;

    dirHandle = CreateFileA(directory,
					FILE_LIST_DIRECTORY, 
					FILE_SHARE_WRITE|FILE_SHARE_READ|FILE_SHARE_DELETE, 
					NULL, 
					OPEN_EXISTING, 
					FILE_FLAG_BACKUP_SEMANTICS|FILE_FLAG_OVERLAPPED, 
                    NULL);

    if (dirHandle != INVALID_HANDLE_VALUE) {
		completionPort = CreateIoCompletionPort(dirHandle, NULL, 0, 0);
		if (completionPort == NULL) {
			CloseHandle(dirHandle);
			dirHandle = INVALID_HANDLE_VALUE;
		}
	}
}

bool FolderWatcher::startListenForChanges() {
	vecFilenames.clear();
	
	if (dirHandle==INVALID_HANDLE_VALUE) return false;
	
	memset(&overlap, 0, sizeof(overlap));
	
	activeBuffer = 1 - activeBuffer;
	
	if (ReadDirectoryChangesW(dirHandle, fileInfoBuffer[activeBuffer], FILE_INFO_BUFFER_SIZE, TRUE, FILE_NOTIFY_CHANGE_LAST_WRITE, NULL, &overlap, NULL)) {
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
	BOOL ok;
	DWORD numBytes, compKey;
	OVERLAPPED *overLapped;

	ok = GetQueuedCompletionStatus(completionPort, &numBytes, &compKey, &overLapped, millis);
	
	if (ok) {
		int oldBuffer = activeBuffer;
		startListenForChanges();
		return processChanges(oldBuffer);
	}
	return 0;
}
	
int FolderWatcher::processChanges(int which) {
	FILE_NOTIFY_INFORMATION *info;

	char *fileInfoBufferPtr = fileInfoBuffer[which];
	
	isListening = false;	
	
    do
    {
		char filename[MAX_PATH];
		info = (FILE_NOTIFY_INFORMATION *) fileInfoBufferPtr;
		
        fileInfoBufferPtr += info->NextEntryOffset;
        int count = WideCharToMultiByte(CP_ACP, 0, info->FileName, info->FileNameLength / sizeof(WCHAR), filename, MAX_PATH - 1, NULL, NULL);
        filename[count] = 0;
		vecFilenames.push_back(filename);
         printf("%i: %s\n", vecFilenames.size(), filename);
	} while (info->NextEntryOffset > 0);
	
	return vecFilenames.size();
}

FolderWatcher::~FolderWatcher() {	
	if (isListening) {
		CancelIo(dirHandle);
	}
	if (completionPort != NULL) CloseHandle(completionPort);
	CloseHandle(dirHandle);
}
