/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
package nl.fcdonders.fieldtrip;

import java.nio.*;

/** A class for wrapping FieldTrip buffer events
*/
public class BufferEvent {
	public BufferEvent() {
		wType  = new WrappedObject();
		wValue = new WrappedObject();
		sample     = -1;
		offset     = 0;
		duration   = 0;
	}
	
	public BufferEvent(String type, String value, int sample) {
		wType  = new WrappedObject(type);
		wValue = new WrappedObject(value);
		this.sample = sample;
		offset = duration = 0;
	}
	
	public BufferEvent(String type, long value, int sample) {
		wType  = new WrappedObject(type);
		wValue = new WrappedObject(value);
		this.sample = sample;
		offset = duration = 0;
	}		
	
	public BufferEvent(String type, int value, int sample) {
		wType  = new WrappedObject(type);
		wValue = new WrappedObject(value);
		this.sample = sample;
		offset = duration = 0;
	}	
	
	public BufferEvent(String type, short value, int sample) {
		wType  = new WrappedObject(type);
		wValue = new WrappedObject(value);
		this.sample = sample;
		offset = duration = 0;
	}	
	
	public BufferEvent(String type, byte value, int sample) {
		wType  = new WrappedObject(type);
		wValue = new WrappedObject(value);
		this.sample = sample;
		offset = duration = 0;
	}	
	
	public BufferEvent(String type, double value, int sample) {
		wType  = new WrappedObject(type);
		wValue = new WrappedObject(value);
		this.sample = sample;
		offset = duration = 0;
	}	

	public BufferEvent(String type, float value, int sample) {
		wType  = new WrappedObject(type);
		wValue = new WrappedObject(value);
		this.sample = sample;
		offset = duration = 0;
	}		
	
	public BufferEvent(ByteBuffer buf) {
		wType  = new WrappedObject();
		wValue = new WrappedObject();	
	
		wType.type   = buf.getInt();
		wType.numel  = buf.getInt();
		wValue.type  = buf.getInt();
		wValue.numel = buf.getInt();
		sample       = buf.getInt();
		offset       = buf.getInt();
		duration     = buf.getInt();
		int size = buf.getInt();
	
		wType.array  = DataType.getObject(wType.type, wType.numel, buf);
		if (wType.array != null) {
			wType.size = wType.numel * DataType.wordSize[wType.type];
		}
	
		wValue.array = DataType.getObject(wValue.type, wValue.numel, buf);
		if (wValue.array != null) {
			wValue.size = wValue.numel * DataType.wordSize[wValue.type];
		}
		
		size -= wType.size + wValue.size;
		if (size != 0) {
			buf.position(buf.position() + size);
		}
	}
	
	public WrappedObject getType() {
		return wType;
	}
	
	public WrappedObject getValue() {
		return wValue;
	}	
	
	public boolean setType(Object typeObj) {
		wType = new WrappedObject(typeObj);
		return wType.type != DataType.UNKNOWN;
	}
	
	public boolean setValue(Object valueObj) {
		wValue = new WrappedObject(valueObj);
		return wValue.type != DataType.UNKNOWN;
	}
	
	public boolean setValueUnsigned(byte[] array) {
		wValue = new WrappedObject();
		wValue.array = array.clone();
		wValue.numel = array.length;
		wValue.size  = array.length;
		wValue.type  = DataType.UINT8;
		return true;
	}

	public int size() {
		return 32 + wType.size + wValue.size;
	}
	
	public static int count(ByteBuffer buf) {
		int num = 0;
		int pos = buf.position();
	
		while (buf.remaining() >= 32) {
			int typeType   = buf.getInt();
			int typeNumEl  = buf.getInt();
			int valueType  = buf.getInt();
			int valueNumEl = buf.getInt();
			buf.getInt(); // sample
			buf.getInt(); // offset
			buf.getInt(); // duration
			int size = buf.getInt();
			int sizeType  = typeNumEl  * DataType.wordSize[typeType];
			int sizeValue = valueNumEl * DataType.wordSize[valueType];
		
			if (sizeType < 0 || sizeValue < 0 || sizeType + sizeValue > size) {
				return -(1+num);
			}
		
			buf.position(buf.position() + size);
			num++;
		}
		buf.position(pos);
		return num;
	}
	
	public void serialize(ByteBuffer buf) {
		buf.putInt(wType.type);
		buf.putInt(wType.numel);
		buf.putInt(wValue.type);
		buf.putInt(wValue.numel);
		buf.putInt(sample);
		buf.putInt(offset);
		buf.putInt(duration);
		buf.putInt(wType.size+wValue.size);
		wType.serialize(buf);
		wValue.serialize(buf);
	}
	
	public void print() {
		System.out.println("type_type   " + wType.type);
		System.out.println("type_numel  " + wType.numel);
		System.out.println("value_type  " + wValue.type);
		System.out.println("value_numel " + wValue.numel);
		System.out.println("size        " + size());
		System.out.println(wType);
		System.out.println(wValue);
	}
	
	public int sample;
	public int offset;
	public int duration;
	protected WrappedObject wType;
	protected WrappedObject wValue;
}