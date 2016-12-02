/** Two small classes for handling configurations for channel selection,
    downsampling, and similar purposes.

    (C) 2010 Stefan Klanke
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include <SignalConfiguration.h>

struct TokenAndSep {
	TokenAndSep() : text(), separator(-1) {};
	TokenAndSep(int len, const char *str, int sep=0) : text(str, len), separator(sep) {};

	bool equals(const char *str) {
		unsigned int i = 0;
		for (i=0;i<text.size();i++) {
			if (toupper(text[i]) != toupper(str[i])) return false;
		}
		if (str[i]!=0) return false;
		return true;
	}

	bool isNumber() {
		if (text.size()==0) return false;
		return isdigit(text[0]);
	}

	std::string text;
	int separator;
};

/** Simple C++ style wrapper around strtol, converts a string to an integer.
	Returns true on success, false on error. The second parameter value receives
	the result of a successful conversion.
*/
bool convertToInt(const std::string& in, int& value) {
	const char *start = in.c_str();
	char *end;
	long v = strtol(start, &end, 10);
	if (start == end) return false;
	if (*end!=0) return false;
	value = (int) v;
	return true;
}

/** Simple C++ style wrapper around strtod, converts a string to a double precision number.
	Returns true on success, false on error. The second parameter value receives
	the result of a successful conversion.
*/
bool convertToDouble(const std::string& in, double& value) {
	const char *start = in.c_str();
	char *end;
	double v = strtod(start, &end);
	if (start == end) return false;
	if (*end!=0) return false;
	value = v;
	return true;
}

bool tokenizeString(int len, const char *str, std::vector<TokenAndSep>& S) {
	TokenAndSep ts;

	int pos = 0;
	int state = 0; // 0=whitespace, 1=in token, 2=in quotes (also in token), 3=find separator
	int startIdx = 0, endIdx = 0; // start and end of the current token

	S.clear();

	for (pos=0;pos<len;pos++) {
		int sp = str[pos];

		if (sp==0 || sp == '\n') break;

		switch (state) {
			case 0: // whitespace
				// continue if this one still is whitespace
				if (isspace(sp)) continue;
				// if token starts with comma or equal sign, this is an error
				if (sp == ',' || sp == '=') return false;
				if (sp == '"') {
					startIdx = pos+1;
					state = 2;
				} else {
					startIdx = pos;
					state = 1;
				}
				break;
			case 1: // token, but not in quotes
				// if there is a quote inside the token, this is an error
				if (sp == '"') return false;
				// comma or equal sign ends the token, have separator
				if (sp==',' || sp=='=') {
					endIdx = pos;
					S.push_back(TokenAndSep(endIdx - startIdx, str + startIdx, sp));
					state  = 0;
				}
				// whitespace ends the token, look for separator
				if (isspace(sp)) {
					endIdx = pos;
					state = 3;
				}
				// else continue in next loop
				break;
			case 2: // token in quotes
				if (sp == '"') {
					// end of token, now look for separator
					endIdx = pos;
					state = 3;
				}
				// else continue in next loop
				break;
			case 3: // look for separator
				if (sp==',' || sp=='=') {
					S.push_back(TokenAndSep(endIdx - startIdx, str + startIdx, sp));
					state  = 0;
				} else if (!isspace(sp)) {
					if (pos > endIdx) {
						S.push_back(TokenAndSep(endIdx - startIdx, str + startIdx, ' '));
						startIdx = pos;
						state  = 1;
					} else {
						// everything but comma,=,white space is an error
						if (!isspace(sp)) return false;
					}
				}
				break;
		}
	}

	// end of loop
	switch(state) {
		case 0:
			return true;
		case 1:
			// if we're inside a token, take pos as endIdx and
			// add the token with \0 as separator
			S.push_back(TokenAndSep(pos - startIdx, str + startIdx, 0));
			return true;
		case 2:
			// if we're still within quotes, this is an error
			return false;
		case 3:
			// if we're looking for a separator, add
			// the token with \0 as separator
			S.push_back(TokenAndSep(endIdx - startIdx, str + startIdx, 0));
			return true;
	}
	return false;
}

bool ChannelSelection::parseString(int len, const char *str) {
	int pos = 0;
	index.clear();
	label.clear();

	while (pos < len) {
		int idx;
		int start;

		// skip white space
		while (isspace(str[pos])) {
			if (++pos == len) return true;
		}

		// next character should be 0-9
		if (!isdigit(str[pos])) return false;
		idx = str[pos] - '0';
		if (++pos == len) return false;

		while (isdigit(str[pos])) {
			idx = 10*idx + (str[pos] - '0');
			if (++pos == len) return false;
		}

		if (idx==0) return false;

		// next character should be =
		if (str[pos] != '=') return false;
		if (++pos == len) return false;

		if (str[pos] == '"') {
			start = ++pos;

			// search next white space
			while (pos < len && str[pos]!='"') pos++;
			// check for empty "" or unterminated "....
			if (pos == len || pos == start) return false;

			// got a label "like this"
			index.push_back(idx-1);
			label.push_back(std::string(str + start, pos-start));
			pos++;
		} else {
			// mark label start
			start = pos;

			// search next white space
			while (!isspace(str[pos]) && pos < len) pos++;
			if (start == pos) return false;

			// got a plain label
			index.push_back(idx-1);
			label.push_back(std::string(str + start, pos-start));
		}
	}
	return true;
}


int SignalConfiguration::parseFile(const char *filename) {
	FILE *fp;
	char line[2048];
	int lineCount = 0;
	int errorCount = 0;
	int addSave = 1, addStream = 1;
    int numTok = 0;

	std::vector<TokenAndSep> TS;

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

		if (!tokenizeString(sizeof(line), line, TS)) goto reportError;

		numTok = TS.size();

		if (numTok == 0) goto reportError;
		if (TS[numTok-1].separator != 0) goto reportError;

		/*
		for (unsigned int i=0;i<TS.size();i++) {
			printf("(%s)%c", TS[i].text.c_str(), TS[i].separator);
		}
		printf("\n");
		*/
		if (numTok == 1) {
			if (TS[0].equals("[select]")) {
				addStream = addSave = 1;
				continue;
			}
			if (TS[0].equals("[stream]")) {
				addStream = 1;
				addSave = 0;
				continue;
			}
			if (TS[0].equals("[save]")) {
				addSave = 1;
				addStream = 0;
				continue;
			}
			goto reportError;
			// no other line with one token only is valid
		}
		if (numTok == 2) {
			if (TS[0].isNumber()) {
				if (TS[0].separator != '=') goto reportError;
				int chn;

				if (!convertToInt(TS[0].text, chn) || chn <= 0) goto reportError;
				if (addSave) {
					if (chn>maxChanSave) maxChanSave = chn;
					chanSelSave.add(chn-1, TS[1].text);
				}
				if (addStream) {
					if (chn>maxChanStream) maxChanStream = chn;
					chanSelStream.add(chn-1, TS[1].text);
				}
				continue;
			}
			if (TS[0].equals("downsample")) {
				int ds;
				if (!convertToInt(TS[1].text, ds) || ds < 1) goto reportError;
				downSample=ds;
				continue;
			}
			if (TS[0].equals("bworder")) {
				int bwOrder;
				if (!convertToInt(TS[1].text, bwOrder) || bwOrder<0) goto reportError;
				order = bwOrder;
				continue;
			}
			if (TS[0].equals("bandwidth")) {
				double bw;
				if (!convertToDouble(TS[1].text, bw) || bw < 0) goto reportError;
				bandwidth = bw;
				continue;
			}
			if (TS[0].equals("samplerate")) {
				double sr;
				// sampleRate = 0 would be interpreted as "leave at default / maximum"
				if (!convertToDouble(TS[1].text, sr) || sr < 0) goto reportError;
				sampleRate = sr;
				continue;
			}
			if (TS[0].equals("batteryrefresh")) {
				int br;
				if (!convertToInt(TS[1].text, br) || br<0) goto reportError;
				batteryRefresh = br;
				continue;
			}
			if (TS[0].equals("statusrefresh")) {
				int sr;
				if (!convertToInt(TS[1].text, sr) || sr<0) goto reportError;
				statusRefresh = sr;
				continue;
			}
			// no other 2-token lines are valid
			goto reportError;
		}
		if (numTok == 3) {
			if (TS[0].equals("splittrigger") && TS[1].separator==',') {
				splitTrigger = true;
				lowTriggerName  = TS[1].text;
				highTriggerName = TS[2].text;
				continue;
			}
			// no other 3-token lines are valid
			goto reportError;
		}
		reportError:

		fprintf(stderr, "Could not parse line %i:\n%s\n", lineCount, line);
		errorCount++;
	}
	fclose(fp);
	return errorCount;
}
