/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
#ifndef __FolderWatcher_h
#define __FolderWatcher_h

#include <vector>
#include <string>
#include <windows.h>

#define FILE_INFO_BUFFER_SIZE  20000

class FolderWatcher {
	public:
	
	FolderWatcher(const char *directory);
	~FolderWatcher();
		
	/* Asynchronous operation, returns "false" on error */
	bool startListenForChanges();
	bool stopListenForChanges();
	
	/* Asynchronous operation, returns number of changes -- if non-zero, restart listening */
	int checkHasChanged(unsigned int ms);
	
	/* retrieve vector of filenames */
	const std::vector<std::string>& getFilenames() const { 
		return vecFilenames; 
	}
	
	bool isValid() const { 
		return (dirHandle != INVALID_HANDLE_VALUE); 
	}
	
	protected:
	
	bool isListening;
	int processChanges(int which);
	
	std::vector<std::string> vecFilenames;
	char fileInfoBuffer[2][FILE_INFO_BUFFER_SIZE];
	int activeBuffer;
	HANDLE dirHandle;
	HANDLE completionPort;
	OVERLAPPED overlap;
};

#endif
