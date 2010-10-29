/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
import java.nio.*;

/** A class for wrapping FieldTrip buffer events
	TODO: provide more constructor to create events from often used
	type/value combinations, provide function for serialization.
*/
class Event {
	public Event(ByteBuffer buf) {
		typeType   = buf.getInt();
		typeNumEl  = buf.getInt();
		valueType  = buf.getInt();
		valueNumEl = buf.getInt();
		sample     = buf.getInt();
		offset     = buf.getInt();
		duration   = buf.getInt();
		int size = buf.getInt();
		int sizeType  = valueNumEl * Header.wordsize[valueType];
		int sizeValue = typeNumEl * Header.wordsize[typeType];
	
		value = new byte[sizeType];
		type  = new byte[sizeValue];
	
		buf.get(type);
		buf.get(value);
		size -= sizeType + sizeValue;
		if (size != 0) {
			buf.position(buf.position() + size);
		}
	}
	
	public String getTypeString() {
		if (typeType != Header.DATATYPE_CHAR) return null;
		return new String(type);
	}
	
	public String getValueString() {
		if (valueType != Header.DATATYPE_CHAR) return null;
		return new String(value);
	}	
	
	public int sample;
	public int offset;
	public int duration;
	public int typeType;
	public int typeNumEl;
	public int valueType;
	public int valueNumEl;
	byte[] value;
	byte[] type;
}