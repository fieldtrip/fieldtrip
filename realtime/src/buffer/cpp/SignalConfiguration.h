/** Two small classes for handling configurations for channel selection,
    downsampling, and similar purposes.

    (C) 2010 Stefan Klanke
 */

#ifndef __SignalConfiguration_h
#define __SignalConfiguration_h

#include <vector>
#include <string>

/** Simple C++ style wrapper around strtol, converts a string to an integer.
	Returns true on success, false on error. The second parameter value receives
	the result of a successful conversion.
*/
bool convertToInt(const std::string& in, int& value);

/** Simple C++ style wrapper around strtod, converts a string to a double precision number.
	Returns true on success, false on error. The second parameter value receives
	the result of a successful conversion.
*/
bool convertToDouble(const std::string& in, double& value);

struct ChannelSelection {
	std::vector<int> index;
	std::vector<std::string> label;

	bool parseString(int len, const char *str);

	void clear() {
		index.clear();
		label.clear();
	}

	void add(int newIndex, const char *newLabel) {
		index.push_back(newIndex);
		label.push_back(newLabel);
	}

	void add(int newIndex, const std::string& newLabel) {
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

	static ChannelSelection selectByLabel(const ChannelSelection& sel, const ChannelSelection& all) {
		ChannelSelection cs;
		for (unsigned int i=0;i<sel.index.size();i++) {
			for (unsigned int j=0;j<all.index.size();j++) {
				if (sel.label[i].compare(all.label[j]) == 0) {
					cs.index.push_back(all.index[j]);
					cs.label.push_back(all.label[j]);
					break;
				}
			}
		}
		return cs;
	}

	int getMinIndex() const {
		if (index.empty()) return -1;
		int idx = index[0];
		for (unsigned int i=1;i<index.size();i++) {
			if (index[i] < idx) idx = index[i];
		}
		return idx;
	}

	int getMaxIndex() const {
		int idx = -1;
		for (unsigned int i=0;i<index.size();i++) {
			if (index[i] > idx) idx = index[i];
		}
		return idx;
	}
};

class SignalConfiguration {
	public:

	SignalConfiguration() : splitTrigger(0), chanSelSave(), chanSelStream(),
							downSample(1), maxChanSave(0), maxChanStream(0),
							order(0), bandwidth(-1.0), sampleRate(0.0),
							batteryRefresh(10), statusRefresh(2) {};
	~SignalConfiguration() {};

	/** Parse given configuration file for options, returns number of errors */
	int parseFile(const char *filename);
	bool selectChannels(const char *str);
	bool setDownsampling(int factor) {
		if (factor < 1) return false;
		downSample = factor;
		return true;
	}

	void setBandwidth(double bandwidth) {
		this->bandwidth=bandwidth;
	}

	void setSampleRate(double sampleRate) {
		this->sampleRate = sampleRate;
	}

	void setOrder(int order) {
		if (order>=0) {
			this->order = order;
		}
	}

	int getDownsampling() const { return downSample; }
	double getBandwidth() const { return bandwidth; }
	double getSampleRate() const { return sampleRate; }
	int getOrder() const { return order; }
	int getBatteryRefresh() const { return batteryRefresh; }
	int getStatusRefresh() const { return statusRefresh; }
	int getMaxSavingChannel() const { return maxChanSave; }
	int getMaxStreamingChannel() const { return maxChanStream; }
	bool useSplittedTrigger() const { return splitTrigger; }
	const char *getLowTriggerName() const { return lowTriggerName.c_str(); }
	const char *getHighTriggerName() const { return highTriggerName.c_str(); }

	const ChannelSelection& getSavingSelection() const { return chanSelSave; }
	const ChannelSelection& getStreamingSelection() const { return chanSelStream; }

	bool setSavingSelection(const ChannelSelection& sel) {
		if (sel.index.size() != sel.label.size()) return false;
		chanSelSave = sel;
		maxChanSave = sel.getMaxIndex();
		return true;
	}

	bool setStreamingSelection(const ChannelSelection& sel) {
		if (sel.index.size() != sel.label.size()) return false;
		chanSelStream = sel;
		maxChanStream = sel.getMaxIndex();
		return true;
	}

	void selectForStreaming(int index, const char *label) {
		chanSelStream.add(index, label);
		if (index > maxChanStream) maxChanStream = index;
	}

	void selectForSaving(int index, const char *label) {
		chanSelSave.add(index, label);
		if (index > maxChanSave) maxChanSave = index;
	}

	protected:

	std::string lowTriggerName, highTriggerName;
	bool splitTrigger;
	ChannelSelection chanSelSave;
	ChannelSelection chanSelStream;
	int downSample;
	int maxChanSave;
	int maxChanStream;
	int order;
	double bandwidth;
	double sampleRate;
	int batteryRefresh, statusRefresh;
};

#endif
