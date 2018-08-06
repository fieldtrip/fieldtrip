package nl.fcdonders.fieldtrip.bufferserver.data;

import java.nio.ByteOrder;
import java.nio.ByteBuffer;

/**
 * Wrapper for passing data between the NetworkProtocol and dataStore.
 *
 * @author Wieke Kanters
 *
 */
public class Data {
	public final byte[][][] data;
	public final int dataType;
	public final int nChans;
	public final int nSamples;
	public final ByteOrder order;

	/**
	 * Constructor.
	 *
	 * @param nChans
	 *            number of channels
	 * @param nSamples
	 *            number of samples
	 * @param nBytes
	 *            number of bytes per data point
	 * @param dataType
	 *            dataType
	 * @param data
	 *            3d (nSamples by nChans by nBytes) byte array containing data
	 * @param order
	 *            endianess of the data.
	 */
	public Data(int nChans, int nSamples, int dataType, byte[][][] data,
			ByteOrder order) {
		this.data = data;
		this.dataType = dataType;
		this.nChans = nChans;
		this.nSamples = nSamples;
		this.order = order;
	}

	/**
	 * Returns the size in bytes.
	 *
	 * @return
	 */
	public int size() {
		return nSamples * nChans;
	}
}
