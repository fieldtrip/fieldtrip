/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
package nl.fcdonders.fieldtrip;

import java.nio.channels.*;
import java.nio.*;
import java.net.*;
import java.io.*;

public class BufferClient {
	public static final short VERSION = 1;
	public static final short GET_HDR = 0x201;
	public static final short GET_DAT = 0x202;
	public static final short GET_EVT = 0x203;
	public static final short GET_OK  = 0x204;
	public static final short GET_ERR = 0x205;

	public static final short PUT_HDR = 0x101;
	public static final short PUT_DAT = 0x102;
	public static final short PUT_EVT = 0x103;
	public static final short PUT_OK  = 0x104;
	public static final short PUT_ERR = 0x105;

	public static final short FLUSH_HDR = 0x301;
	public static final short FLUSH_DAT = 0x302;
	public static final short FLUSH_EVT = 0x303;
	public static final short FLUSH_OK  = 0x304;
	public static final short FLUSH_ERR = 0x305;

	public static final short WAIT_DAT = 0x402;
	public static final short WAIT_OK  = 0x404;
	public static final short WAIT_ERR = 0x405;


	public BufferClient() {
		myOrder = ByteOrder.nativeOrder();
	}
	
	public BufferClient(ByteOrder order) {
		myOrder = order;
	}
	
	public boolean connect(String hostname, int port) throws IOException {
		if (sockChan == null) {
			sockChan = SocketChannel.open();		
		} else if (sockChan.isConnected()) {
			disconnect();
		}
		sockChan.connect(new InetSocketAddress(hostname, port));
	
		return sockChan.isConnected();
	}
	
	public boolean connect(String address) throws IOException {
		int colonPos = address.lastIndexOf(':');
		if (colonPos != -1) {
			String hostname = address.substring(0,colonPos);
			Integer port;
		
			try {
				port = new Integer(address.substring(colonPos+1));
			}
			catch (NumberFormatException e) {
				System.out.println(e);
				return false;
			}
			return connect(hostname, port.intValue());
		}
		System.out.println("Address format not recognized / supported yet.");
		// other addresses not recognised yet
		return false;
	}
	
	public void disconnect() throws IOException {
		sockChan.socket().close();
		sockChan = null;
	}
	
	public boolean isConnected() {
		if (sockChan == null) return false;
		return sockChan.isConnected();
	}
	
	public Header getHeader() throws IOException {
		ByteBuffer buf;

		buf = ByteBuffer.allocate(8);
		buf.order(myOrder);
	
		buf.putShort(VERSION).putShort(GET_HDR).putInt(0).rewind();
		writeAll(buf);
	
		buf = readResponse(GET_OK);
		return new Header(buf);
	}
	
	/** Returns true if channel names were written */
	public boolean putHeader(Header hdr) throws IOException {
		ByteBuffer buf;
		int bufsize = hdr.getSerialSize();

		buf = ByteBuffer.allocate(8 + bufsize);
		buf.order(myOrder);
	
		buf.putShort(VERSION).putShort(PUT_HDR).putInt(bufsize);
		hdr.serialize(buf);
		buf.rewind();
		writeAll(buf);
		readResponse(PUT_OK);
		return hdr.channelNameSize > hdr.nChans;
	}
	
	public short[][] getShortData(int first, int last) throws IOException {
		DataDescription dd = new DataDescription();
		ByteBuffer buf = getRawData(first, last, dd);
	
		int nSamples = dd.nSamples;
		int nChans = dd.nChans;
		
		short[][] data = new short[nSamples][nChans];
		
		switch (dd.dataType) {
			case DataType.INT8:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (short) buf.get();
					}
				}
				break;
			case DataType.INT16:
				ShortBuffer sBuf = buf.asShortBuffer();
				for (int n=0;n<nSamples;n++) sBuf.get(data[n]);
				break;
			default:
				System.out.println("Not supported yet - returning zeros.");
		}
	
		return data;
	}
	
	public int[][] getIntData(int first, int last) throws IOException {
		DataDescription dd = new DataDescription();
		ByteBuffer buf = getRawData(first, last, dd);
	
		int nSamples = dd.nSamples;
		int nChans = dd.nChans;
		
		int[][] data = new int[nSamples][nChans];
		
		switch (dd.dataType) {
			case DataType.INT8:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (int) buf.get();
					}
				}
				break;
			case DataType.INT16:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (int) buf.getShort();
					}
				}
				break;
			case DataType.INT32:
				IntBuffer iBuf = buf.asIntBuffer();
				for (int n=0;n<nSamples;n++) iBuf.get(data[n]);
				break;
			default:
				System.out.println("Not supported yet - returning zeros.");
		}
	
		return data;
	}
	
	public long[][] getLongData(int first, int last) throws IOException {
		DataDescription dd = new DataDescription();
		ByteBuffer buf = getRawData(first, last, dd);
	
		int nSamples = dd.nSamples;
		int nChans = dd.nChans;
		
		long[][] data = new long[nSamples][nChans];
		
		switch (dd.dataType) {
			case DataType.INT8:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (int) buf.get();
					}
				}
				break;
			case DataType.INT16:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (int) buf.getShort();
					}
				}
				break;
			case DataType.INT32:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (int) buf.getInt();
					}
				}
				break;
			case DataType.INT64:
				LongBuffer lBuf = buf.asLongBuffer();
				for (int n=0;n<nSamples;n++) lBuf.get(data[n]);
				break;
			default:
				System.out.println("Not supported yet - returning zeros.");
		}
	
		return data;
	}		
	
	public float[][] getFloatData(int first, int last) throws IOException {
		DataDescription dd = new DataDescription();
		ByteBuffer buf = getRawData(first, last, dd);
	
		int nSamples = dd.nSamples;
		int nChans = dd.nChans;
		
		float[][] data = new float[nSamples][nChans];
		
		switch (dd.dataType) {
			case DataType.INT8:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (float) buf.get();
					}
				}
				break;
			case DataType.INT16:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (float) buf.getShort();
					}
				}
				break;
			case DataType.INT32:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (float) buf.getInt();
					}
				}
				break;
			case DataType.FLOAT32:
				FloatBuffer fBuf = buf.asFloatBuffer();
				for (int n=0;n<nSamples;n++) fBuf.get(data[n]);
				break;
			case DataType.FLOAT64:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (float) buf.getDouble();
					}
				}
				break;
			default:
				System.out.println("Not supported yet - returning zeros.");
		}
	
		return data;
	}
	
	public double[][] getDoubleData(int first, int last) throws IOException {
		DataDescription dd = new DataDescription();
		ByteBuffer buf = getRawData(first, last, dd);
	
		int nSamples = dd.nSamples;
		int nChans = dd.nChans;
		
		double[][] data = new double[nSamples][nChans];
		
		switch (dd.dataType) {
			case DataType.INT8:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (double) buf.get();
					}
				}
				break;
			case DataType.INT16:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (double) buf.getShort();
					}
				}
				break;
			case DataType.INT32:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (double) buf.getInt();
					}
				}
				break;		
			case DataType.INT64:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = (double) buf.getLong();
					}
				}
				break;		
			case DataType.FLOAT32:
				for (int i=0;i<nSamples;i++) {
					for (int j=0;j<nChans;j++) {
						data[i][j] = buf.getFloat();
					}
				}
				break;
			case DataType.FLOAT64:
				DoubleBuffer dBuf = buf.asDoubleBuffer();
				for (int n=0;n<nSamples;n++) dBuf.get(data[n]);
				break;
			default:
				System.out.println("Not supported yet - returning zeros.");
		}
	
		return data;
	}
	
	
	public ByteBuffer getRawData(int first, int last, DataDescription descr) throws IOException {
		ByteBuffer buf;

		buf = ByteBuffer.allocate(16);
		buf.order(myOrder);
	
		buf.putShort(VERSION).putShort(GET_DAT).putInt(8);
		buf.putInt(first).putInt(last).rewind();
		writeAll(buf);
		buf = readResponse(GET_OK);
		
		descr.nChans    = buf.getInt();
		descr.nSamples  = buf.getInt();
		descr.dataType  = buf.getInt();
		descr.sizeBytes = buf.getInt();
	
		int dataSize = descr.nChans * descr.nSamples * DataType.wordSize[descr.dataType];
	
		if (dataSize > descr.sizeBytes || descr.sizeBytes > buf.remaining()) {
			throw new IOException("Invalid size definitions in response from GET DATA request");
		}
	
		return buf.slice();
	}	
	
	
	public BufferEvent[] getEvents() throws IOException {
		ByteBuffer buf;

		buf = ByteBuffer.allocate(8);
		buf.order(myOrder); 
	
		buf.putShort(VERSION).putShort(GET_EVT).putInt(0).rewind();
	
		writeAll(buf);
		buf = readResponse(GET_OK);
	
		int numEvt = BufferEvent.count(buf);
		if (numEvt < 0) throw new IOException("Invalid event definitions in response.");
	
		BufferEvent[] evs = new BufferEvent[numEvt];
		for (int n=0;n<numEvt;n++) {
			evs[n] = new BufferEvent(buf);
		}
		return evs;
	}	
	
	
	public BufferEvent[] getEvents(int first, int last) throws IOException {
		ByteBuffer buf;

		buf = ByteBuffer.allocate(16);
		buf.order(myOrder); 
	
		buf.putShort(VERSION).putShort(GET_EVT).putInt(8);
		buf.putInt(first).putInt(last).rewind();
	
		writeAll(buf);
		buf = readResponse(GET_OK);
	
		int numEvt = BufferEvent.count(buf);
		if (numEvt != (last-first+1)) throw new IOException("Invalid event definitions in response.");
	
		BufferEvent[] evs = new BufferEvent[numEvt];
		for (int n=0;n<numEvt;n++) {
			evs[n] = new BufferEvent(buf);
		}
		return evs;
	}
	
	public void putRawData(int nSamples, int nChans, int dataType, byte[] data) throws IOException {
		if (nSamples == 0) return;
		if (nChans == 0) return;
	
		if (data.length != nSamples*nChans*DataType.wordSize[dataType]) {
			throw new IOException("Raw buffer does not match data description");
		}
		
		ByteBuffer buf = preparePutData(nChans, nSamples, dataType);
		buf.put(data);
		buf.rewind();
		writeAll(buf);
		readResponse(PUT_OK);
	}	
	
	public void putData(byte[][] data) throws IOException {
		int nSamples = data.length;
		if (nSamples == 0) return;
		int nChans = data[0].length;
		if (nChans == 0) return;
	
		for (int i=1;i<nSamples;i++) {
			if (nChans != data[i].length) {
				throw new IOException("Cannot write non-rectangular data array");
			}
		}
		
		ByteBuffer buf = preparePutData(nChans, nSamples, DataType.INT8);
		for (int i=0;i<nSamples;i++) {
			buf.put(data[i]);
		}
		buf.rewind();
		writeAll(buf);
		readResponse(PUT_OK);
	}
	
	public void putData(short[][] data) throws IOException {
		int nSamples = data.length;
		if (nSamples == 0) return;
		int nChans = data[0].length;
		if (nChans == 0) return;
	
		for (int i=1;i<nSamples;i++) {
			if (nChans != data[i].length) {
				throw new IOException("Cannot write non-rectangular data array");
			}
		}
		
		ByteBuffer buf = preparePutData(nChans, nSamples, DataType.INT16);
		for (int i=0;i<nSamples;i++) {
			buf.asShortBuffer().put(data[i]);
		}
		buf.rewind();
		writeAll(buf);
		readResponse(PUT_OK);
	}	
	
	public void putData(int[][] data) throws IOException {
		int nSamples = data.length;
		if (nSamples == 0) return;
		int nChans = data[0].length;
		if (nChans == 0) return;
	
		for (int i=1;i<nSamples;i++) {
			if (nChans != data[i].length) {
				throw new IOException("Cannot write non-rectangular data array");
			}
		}
		
		ByteBuffer buf = preparePutData(nChans, nSamples, DataType.INT32);
		for (int i=0;i<nSamples;i++) {
			buf.asIntBuffer().put(data[i]);
		}
		buf.rewind();
		writeAll(buf);
		readResponse(PUT_OK);
	}
	
	public void putData(long[][] data) throws IOException {
		int nSamples = data.length;
		if (nSamples == 0) return;
		int nChans = data[0].length;
		if (nChans == 0) return;
	
		for (int i=1;i<nSamples;i++) {
			if (nChans != data[i].length) {
				throw new IOException("Cannot write non-rectangular data array");
			}
		}
		
		ByteBuffer buf = preparePutData(nChans, nSamples, DataType.INT64);
		for (int i=0;i<nSamples;i++) {
			buf.asLongBuffer().put(data[i]);
		}
		buf.rewind();
		writeAll(buf);
		readResponse(PUT_OK);
	}	
		
	public void putData(float[][] data) throws IOException {
		int nSamples = data.length;
		if (nSamples == 0) return;
		int nChans = data[0].length;
		if (nChans == 0) return;
	
		for (int i=1;i<nSamples;i++) {
			if (nChans != data[i].length) {
				throw new IOException("Cannot write non-rectangular data array");
			}
		}
		
		ByteBuffer buf = preparePutData(nChans, nSamples, DataType.FLOAT32);
		for (int i=0;i<nSamples;i++) {
			buf.asFloatBuffer().put(data[i]);
		}
		buf.rewind();
		writeAll(buf);
		readResponse(PUT_OK);
	}	
	
	public void putData(double[][] data) throws IOException {
		int nSamples = data.length;
		if (nSamples == 0) return;
		int nChans = data[0].length;
		if (nChans == 0) return;
	
		for (int i=1;i<nSamples;i++) {
			if (nChans != data[i].length) {
				throw new IOException("Cannot write non-rectangular data array");
			}
		}
		
		ByteBuffer buf = preparePutData(nChans, nSamples, DataType.FLOAT64);
		for (int i=0;i<nSamples;i++) {
			buf.asDoubleBuffer().put(data[i]);
		}
		buf.rewind();
		writeAll(buf);
		readResponse(PUT_OK);
	}	

	public void putEvent(BufferEvent e) throws IOException {
		ByteBuffer buf;

		buf = ByteBuffer.allocate(8+e.size());
		buf.order(myOrder); 
	
		buf.putShort(VERSION).putShort(PUT_EVT).putInt(e.size());
		e.serialize(buf);
		buf.rewind();
		writeAll(buf);
		readResponse(PUT_OK);
	}

	public void putEvents(BufferEvent[] e) throws IOException {
		ByteBuffer buf;
		int bufsize = 0;
	
		for (int i=0;i<e.length;i++) {
			bufsize += e[i].size();
		}

		buf = ByteBuffer.allocate(8+bufsize);
		buf.order(myOrder); 
	
		buf.putShort(VERSION).putShort(PUT_EVT).putInt(bufsize);
		for (int i=0;i<e.length;i++) {
			e[i].serialize(buf);
		}
		buf.rewind();
		writeAll(buf);
		readResponse(PUT_OK);
	}		
	
	public void flushHeader() throws IOException {
		ByteBuffer buf = ByteBuffer.allocate(8);
		buf.order(myOrder); 
	
		buf.putShort(VERSION).putShort(FLUSH_HDR).putInt(0).rewind();
		writeAll(buf);
		buf = readResponse(FLUSH_OK);
	}	
	
	public void flushData() throws IOException {
		ByteBuffer buf = ByteBuffer.allocate(8);
		buf.order(myOrder); 
	
		buf.putShort(VERSION).putShort(FLUSH_DAT).putInt(0).rewind();
		writeAll(buf);
		buf = readResponse(FLUSH_OK);
	}	
	
	public void flushEvents() throws IOException {
		ByteBuffer buf = ByteBuffer.allocate(8);
		buf.order(myOrder); 
	
		buf.putShort(VERSION).putShort(FLUSH_EVT).putInt(0).rewind();
		writeAll(buf);
		buf = readResponse(FLUSH_OK);
	}	
	
	public SamplesEventsCount wait(int nSamples, int nEvents, int timeout) throws IOException {
		ByteBuffer buf;

		buf = ByteBuffer.allocate(20);
		buf.order(myOrder); 
	
		buf.putShort(VERSION).putShort(WAIT_DAT).putInt(12);
		buf.putInt(nSamples).putInt(nEvents).putInt(timeout).rewind();
	
		writeAll(buf);
		buf = readResponse(WAIT_OK);
	
		return new SamplesEventsCount(buf.getInt(), buf.getInt());
	}
	
	public SamplesEventsCount waitForSamples(int nSamples, int timeout) throws IOException {
		return wait(nSamples, -1, timeout);
	}	
	
	public SamplesEventsCount waitForEvents(int nEvents, int timeout) throws IOException {
		return wait(-1, nEvents, timeout);
	}		
	
	public SamplesEventsCount poll() throws IOException {
		return wait(0,0,0);
	}
	
	//*********************************************************************
	//		protected methods and variables from here on
	//*********************************************************************
	
	protected ByteBuffer readAll(ByteBuffer dst) throws IOException {
		int rem = dst.remaining();
		while (rem > 0) {
			int now = sockChan.read(dst);
			rem -= now;
		}
		return dst;
	}
	
	protected ByteBuffer readResponse(int expected) throws IOException {
		ByteBuffer def = ByteBuffer.allocate(8);
		def.order(myOrder);
		readAll(def);
		def.rewind();
	
		if (def.getShort() != VERSION) throw new IOException("Invalid VERSION returned.");
		if (def.getShort() != expected) throw new IOException("Error returned from FieldTrip buffer server.");
		int size = def.getInt();
	
		ByteBuffer buf = ByteBuffer.allocate(size);
		buf.order(myOrder);
		readAll(buf);
		buf.rewind();
		return buf;
	}	
	
	protected ByteBuffer writeAll(ByteBuffer dst) throws IOException {
		int rem = dst.remaining();
		while (rem > 0) {
			int now = sockChan.write(dst);
			rem -= now;
		}
		return dst;
	}	
	
	protected ByteBuffer preparePutData(int nChans, int nSamples, int type) {
		int bufsize = DataType.wordSize[type]*nSamples*nChans;
		
		ByteBuffer buf = ByteBuffer.allocate(8+16+bufsize);
		buf.order(myOrder);
		buf.putShort(VERSION).putShort(PUT_DAT).putInt(16+bufsize);
		buf.putInt(nChans).putInt(nSamples).putInt(type).putInt(bufsize);
		return buf;
	}
	
	protected SocketChannel sockChan;
	protected ByteOrder myOrder;
}