/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
import java.nio.*;

/** A class for wrapping a FieldTrip buffer header structure.
	TODO: also handle chunks other than "channel names",
	provide function for serialization (for writing).
*/
class Header {
	public static final int DATATYPE_CHAR   = 0;
	public static final int DATATYPE_UINT8  = 1;
	public static final int DATATYPE_INT8   = 2;
	public static final int DATATYPE_UINT16 = 3;
	public static final int DATATYPE_INT16  = 4;
	public static final int DATATYPE_UINT32 = 5;
	public static final int DATATYPE_INT32  = 6;
	public static final int DATATYPE_UINT64 = 7;
	public static final int DATATYPE_INT64  = 8;
	public static final int DATATYPE_FLOAT32= 9;
	public static final int DATATYPE_FLOAT64= 10;

	public static final int[] wordsize = {1,1,1,2,2,4,4,8,8,4,8};

	public static final int CHUNK_UNKNOWN = 0;
	public static final int CHUNK_CHANNEL_NAMES = 1;
	public static final int CHUNK_CHANNEL_FLAGS = 2;
	public static final int CHUNK_RESOLUTIONS = 3;
	public static final int CHUNK_ASCII_KEYVAL = 4;
	public static final int CHUNK_NIFTI1 = 5;
	public static final int CHUNK_SIEMENS_AP = 6;
	public static final int CHUNK_CTF_RES4 = 7;
	

	public Header(ByteBuffer buf) {
		nChans   = buf.getInt();
		nSamples = buf.getInt();
		nEvents  = buf.getInt();
		fSample  = buf.getFloat();
		dataType = buf.getInt();
		int size = buf.getInt();
		labels   = new String[nChans];
	
		while (size > 0) {
			int chunkType = buf.getInt();
			int chunkSize = buf.getInt();
			byte[] bs = new byte[chunkSize];
			buf.get(bs);
			
			if (chunkType == CHUNK_CHANNEL_NAMES) {
				int n = 0, len = 0;
				for (int pos = 0;pos<chunkSize;pos++) {
					if (bs[pos]==0) {
						if (len>0) {
							labels[n] = new String(bs, pos-len, len);
						}
						len = 0;
						if (++n == nChans) break;
					} else {
						len++;
					}
				}
			} else {
				// ignore all other chunks for now
			}
			size -= 8 + chunkSize;
		}
	}
	
	public Header(int nChans, float fSample, int dataType) {
		this.nChans   = nChans;
		this.fSample  = fSample;
		this.nSamples = 0;
		this.nEvents  = 0;
		this.dataType = dataType;
		this.labels   = new String[nChans]; // allocate, but do not fill
	}
	
	public int dataType;
	public float fSample;
	public int nChans;
	public int nSamples;
	public int nEvents;
	public String[] labels;
}