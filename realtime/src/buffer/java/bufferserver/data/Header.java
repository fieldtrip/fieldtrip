package nl.fcdonders.fieldtrip.bufferserver.data;

import java.nio.ByteOrder;
import java.nio.ByteBuffer;

/**
 * Wrapper class for header information.
 * 
 * @author lkanters
 * 
 */
public class Header {
	public final int dataType;
	public final float fSample;
	public final int nChans;
	public final int nSamples;
	public final int nEvents;
	public final Chunk[] chunks;
	public final int nChunks;
	public final ByteOrder order;

	/**
	 * Constructor
	 * 
	 * @param header
	 *            A header
	 * @param chunks
	 *            An array of extended header chunks
	 * @param order
	 *            The ByteOrder of the header
	 */
	public Header(Header header, Chunk[] chunks, ByteOrder order) {
		nChans = header.nChans;
		fSample = header.fSample;
		nSamples = header.nSamples;
		nEvents = header.nEvents;
		dataType = header.dataType;
		this.chunks = chunks;
		nChunks = chunks.length;
		this.order = order;
	}

	/**
	 * Constructor
	 * 
	 * @param header
	 *            A header
	 * @param nSamples
	 *            Number of samples
	 * @param nEvents
	 *            Number of events
	 */
	public Header(Header header, int nSamples, int nEvents) {
		nChans = header.nChans;
		fSample = header.fSample;
		this.nSamples = nSamples;
		this.nEvents = nEvents;
		dataType = header.dataType;
		chunks = header.chunks;
		nChunks = header.nChunks;
		order = header.order;
	}

	/**
	 * Constructor
	 * 
	 * @param nChans
	 *            Number of channels
	 * @param fSample
	 *            Sampling frequency
	 * @param dataType
	 *            Datatype of the data
	 * @param order
	 *            Bytorder of the data
	 */
	public Header(int nChans, float fSample, int dataType, ByteOrder order) {
		this.nChans = nChans;
		this.fSample = fSample;
		nSamples = 0;
		nEvents = 0;
		this.dataType = dataType;
		chunks = null;
		nChunks = 0;
		this.order = order;
	}

	/**
	 * Constructor
	 * 
	 * @param nChans
	 *            Number of channels
	 * @param fSample
	 *            Sampling frequency
	 * @param dataType
	 *            Datatype of the data
	 * @param chunks
	 *            Extended header chunks
	 * @param order
	 *            Bytorder of the data
	 */
	public Header(int nChans, float fSample, int dataType, Chunk[] chunks,
			ByteOrder order) {
		this.nChans = nChans;
		this.fSample = fSample;
		nSamples = 0;
		nEvents = 0;
		this.dataType = dataType;
		this.chunks = chunks;
		nChunks = chunks.length;
		this.order = order;
	}

	/**
	 * Constructor
	 * 
	 * @param nChans
	 *            Number of channels
	 * @param nSamples
	 *            Number of samples
	 * @param nEvents
	 *            Number of events
	 * @param fSample
	 *            Sampling frequency
	 * @param dataType
	 *            Datatype of the data
	 * @param order
	 *            Bytorder of the data
	 */
	public Header(int nChans, int nSamples, int nEvents, int fSample,
			int dataType, ByteOrder order) {
		this.nChans = nChans;
		this.fSample = fSample;
		this.nSamples = nSamples;
		this.nEvents = nEvents;
		this.dataType = dataType;
		chunks = null;
		nChunks = 0;
		this.order = order;
	}

	/**
	 * Constructor
	 * 
	 * @param nChans
	 *            Number of channels
	 * @param nSamples
	 *            Number of samples
	 * @param nEvents
	 *            Number of events
	 * @param fSample
	 *            Sampling frequency
	 * @param dataType
	 *            Datatype of the data
	 * @param chunks
	 *            Extended header chunks
	 * @param order
	 *            Bytorder of the data
	 */
	public Header(int nChans, int nSamples, int nEvents, int fSample,
			int dataType, Chunk[] chunks, ByteOrder order) {
		this.nChans = nChans;
		this.fSample = fSample;
		this.nSamples = nSamples;
		this.nEvents = nEvents;
		this.dataType = dataType;
		this.chunks = chunks;
		nChunks = chunks.length;
		this.order = order;
	}
	 
	protected void serialize(ByteBuffer buf) {
		buf.putInt(nChans);
		buf.putInt(nSamples);
		buf.putInt(nEvents);
		buf.putFloat(fSample);
		buf.putInt(dataType);
		if( chunks.length<=0 ) {
			 buf.putInt(0);
		} else {
			 for ( int i=0; i < chunks.length; i++){
				  buf.putInt(8+chunks[i].size);
				  chunks[i].serialize(buf);
			 }		
		}
	}
}
