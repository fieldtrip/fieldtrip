#ifndef __SignalConfiguration_h
#define __SignalConfiguration_h

#include <vector>
#include <string>

struct ChannelSelection {
	std::vector<int> index;
	std::vector<std::string> label;
	
	void clear() {
		index.clear();
		label.clear();
	}
	void add(int newIndex, const char *newLabel) {
		index.push_back(newIndex);
		label.push_back(newLabel);
	}
	int getSize() const {
		return index.size();
	}
	int getIndex(int n) const {
		return index[n];
	}
	const char *getLabel(unsigned int n) const {
		return label[n].c_str();
	}	
};

class SignalConfiguration {
	public:
	
	SignalConfiguration() : chanSelSave(), chanSelStream(), downSample(1), maxChanSave(0), maxChanStream(0) {};
	~SignalConfiguration() {};
	
	/** Parse given configuration file for options, returns number of errors */
	int parseFile(const char *filename);
	bool selectChannels(const char *str);
	bool setDownsampling(int factor) {
		if (factor < 1) return false;
		downSample = factor;
		return true;
	}
	
	int getDownsampling() const { return downSample; };
	int getMaxSavingChannel() const { return maxChanSave; };
	int getMaxStreamingChannel() const { return maxChanStream; };
	
	const ChannelSelection& getSavingSelection() const { return chanSelSave; };
	const ChannelSelection& getStreamingSelection() const { return chanSelStream; };
	
	protected:
	
	ChannelSelection chanSelSave;
	ChannelSelection chanSelStream;
	int downSample;
	int maxChanSave;
	int maxChanStream;
};

#endif
