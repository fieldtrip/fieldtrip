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


	public BufferClient() throws IOException {
		myOrder = ByteOrder.nativeOrder();
	}
	
	public BufferClient(ByteOrder order) throws IOException {
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
	
	public short[][] getShortData(int first, int last) throws IOException {
		ByteBuffer buf;

		buf = ByteBuffer.allocate(16);
		buf.order(myOrder);
	
		buf.putShort(VERSION).putShort(GET_DAT).putInt(8);
		buf.putInt(first).putInt(last).rewind();
		writeAll(buf);
		buf = readResponse(GET_OK);
	
		int nChans = buf.getInt();
		int nSamples = buf.getInt();
		int dataType = buf.getInt();
		int datSize = buf.getInt();
		
		short[][] data = new short[nSamples][nChans];
		
		switch (dataType) {
			case DataType.UINT8:
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
		ByteBuffer buf;

		buf = ByteBuffer.allocate(16);
		buf.order(myOrder);
	
		buf.putShort(VERSION).putShort(GET_DAT).putInt(8);
		buf.putInt(first).putInt(last).rewind();
		writeAll(buf);
		buf = readResponse(GET_OK);
	
		int nChans = buf.getInt();
		int nSamples = buf.getInt();
		int dataType = buf.getInt();
		int datSize = buf.getInt();
		
		int[][] data = new int[nSamples][nChans];
		
		switch (dataType) {
			case DataType.UINT8:
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
	
	public float[][] getFloatData(int first, int last) throws IOException {
		ByteBuffer buf;

		buf = ByteBuffer.allocate(16);
		buf.order(myOrder);
	
		buf.putShort(VERSION).putShort(GET_DAT).putInt(8);
		buf.putInt(first).putInt(last).rewind();
		writeAll(buf);
		buf = readResponse(GET_OK);
	
		int nChans = buf.getInt();
		int nSamples = buf.getInt();
		int dataType = buf.getInt();
		int datSize = buf.getInt();
		
		float[][] data = new float[nSamples][nChans];
		
		switch (dataType) {
			case DataType.UINT8:
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
		ByteBuffer buf;

		buf = ByteBuffer.allocate(16);
		buf.order(myOrder);
	
		buf.putShort(VERSION).putShort(GET_DAT).putInt(8);
		buf.putInt(first).putInt(last).rewind();
		writeAll(buf);
		buf = readResponse(GET_OK);
	
		int nChans = buf.getInt();
		int nSamples = buf.getInt();
		int dataType = buf.getInt();
		int datSize = buf.getInt();
		
		double[][] data = new double[nSamples][nChans];
		
		switch (dataType) {
			case DataType.UINT8:
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
	
	protected SocketChannel sockChan;
	protected ByteOrder myOrder;
}