/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
package nl.fcdonders.fieldtrip;

import java.nio.*;

/** A class for wrapping FieldTrip buffer events
	TODO: provide more constructor to create events from often used
	type/value combinations, provide function for serialization.
*/
public class BufferEvent {
	public BufferEvent(ByteBuffer buf) {
		typeType   = buf.getInt();
		typeNumEl  = buf.getInt();
		valueType  = buf.getInt();
		valueNumEl = buf.getInt();
		sample     = buf.getInt();
		offset     = buf.getInt();
		duration   = buf.getInt();
		int size = buf.getInt();
		int sizeType  = typeNumEl * DataType.wordSize[typeType];
		int sizeValue = valueNumEl * DataType.wordSize[valueType];	
	
		type  = DataType.getObject(typeType, typeNumEl, buf);
		value = DataType.getObject(valueType, valueNumEl, buf);
	
		size -= sizeType + sizeValue;
		if (size != 0) {
			buf.position(buf.position() + size);
		}
	}
	
	public Object getType() {
		return type;
	}
	
	public Object getValue() {
		return value;
	}		
		
	public void setType(Object typeObj) {
		type = typeObj;
	}
	
	public void setValue(Object valueObj) {
		value = valueObj;
	}
		
	public int size() {
		int sizeType  = typeNumEl * DataType.wordSize[typeType];
		int sizeValue = valueNumEl * DataType.wordSize[valueType];
	
		return 32 + sizeType + sizeValue;
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
			int sizeType  = typeNumEl * DataType.wordSize[typeType];
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
	
	
	public int sample;
	public int offset;
	public int duration;
	public int typeType;
	public int typeNumEl;
	public int valueType;
	public int valueNumEl;
	Object value;
	Object type;
}