/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
package nl.fcdonders.fieldtrip;

import java.nio.*;

/** A class for wrapping relevant Java objects in a way that
	is easily convertible to FieldTrip data types.
*/
public class WrappedObject {
	protected int type;
	protected int numel;
	protected int size;
	protected Object array;
	
	public WrappedObject() {
		type  = DataType.UNKNOWN;
		numel = 0;
		size  = 0;
		array = null;
	}
	
	public WrappedObject(String s) {
		type  = DataType.CHAR;
		numel = s.getBytes().length;
		size  = numel;
		array = s;
	}
	
	public WrappedObject(double x) {
		type  = DataType.FLOAT64;
		numel = 1;
		size  = 8;
		array = new double[] {x};
	}	
	
	public WrappedObject(float x) {
		type  = DataType.FLOAT32;
		numel = 1;
		size  = 4;
		array = new float[] {x};
	}
	
	public WrappedObject(long x) {
		type  = DataType.INT64;
		numel = 1;
		size  = 8;
		array = new long[] {x};
	}	
	
	public WrappedObject(int x) {
		type  = DataType.INT32;
		numel = 1;
		size  = 4;
		array = new int[] {x};
	}
	
	public WrappedObject(short x) {
		type  = DataType.INT16;
		numel = 1;
		size  = 2;
		array = new short[] {x};
	}
	
	public WrappedObject(byte x) {
		type  = DataType.INT8;
		numel = 1;
		size  = 1;
		array = new byte[] {x};
	}	
	
	public WrappedObject(Object obj) {
		this();
	
		Class cls = obj.getClass();
		String name = cls.getName();

		if (cls.isArray()) {
			Class elc = cls.getComponentType();
			if (!elc.isPrimitive()) return;
			
			if        (name == "[D") {
				type = DataType.FLOAT64;
				array = ((double[]) obj).clone();
				numel = ((double[]) obj).length;
			} else if (name == "[F") {
				type = DataType.FLOAT32;
				array = ((float[]) obj).clone();
				numel = ((float[]) obj).length;
			} else if (name == "[J") {
				type = DataType.INT64;
				array = ((long[]) obj).clone();
				numel = ((long[]) obj).length;
			} else if (name == "[I") {
				type = DataType.INT32;
				array = ((int[]) obj).clone();
				numel = ((int[]) obj).length;
			
			} else if (name == "[S") {
				type = DataType.INT16;
				array = ((short[]) obj).clone();
				numel = ((short[]) obj).length;
			
			} else if (name == "[B") {
				type = DataType.INT8;
				array = ((byte[]) obj).clone();
				numel = ((byte[]) obj).length;
			
			} else {
				return; // keep as unknown
			}
			size  = numel * DataType.wordSize[type];
			return;
		} else if (name == "java.lang.String") {
			type = DataType.CHAR;
			array = obj;
			numel = ((String) obj).getBytes().length;
			size  = numel;
			return;
		} else if (name == "java.lang.Double") {
			type = DataType.FLOAT64;
			array = new double[] {((Double) obj).doubleValue()};
		} else if (name == "java.lang.Float") {
			type = DataType.FLOAT32;
			array = new float[] {((Float) obj).floatValue()};
		} else if (name == "java.lang.Long") {
			type = DataType.INT64;
			array = new long[] {((Long) obj).longValue()};
		} else if (name == "java.lang.Integer") {
			type = DataType.INT32;
			array = new int[] {((Integer) obj).intValue()};
		} else if (name == "java.lang.Short") {
			type = DataType.INT16;
			array = new short[] {((Short) obj).shortValue()};
		} else if (name == "java.lang.Byte") {
			type = DataType.INT8;
			array = new byte[] {((Byte) obj).byteValue()};		
		} else {
			return;
		}
		numel = 1;
		size  = DataType.wordSize[type];
	}	
		
	public void serialize(ByteBuffer buf) {
		switch(type) {
			case DataType.CHAR:
				buf.put(((String) array).getBytes());
				break;
			case DataType.UINT8:
			case DataType.INT8:
				buf.put((byte[]) array);
				break;
			case DataType.UINT16:
			case DataType.INT16:
				buf.asShortBuffer().put((short[]) array);
				break;
			case DataType.UINT32:
			case DataType.INT32:
				buf.asIntBuffer().put((int[]) array);
				break;
			case DataType.UINT64:
			case DataType.INT64:
				buf.asLongBuffer().put((long[]) array);
				break;
			case DataType.FLOAT32:
				buf.asFloatBuffer().put((float[]) array);
				break;
			case DataType.FLOAT64:
				buf.asDoubleBuffer().put((double[]) array);
				break;
		}
	}	
	
	public String toString() {
		if (type == DataType.CHAR) return (String) array;
		if (type == DataType.FLOAT64) return String.valueOf(((double[]) array)[0]);
		if (type == DataType.FLOAT32) return String.valueOf(((float[]) array)[0]);
		if (type == DataType.INT64) return String.valueOf(((long[]) array)[0]);
		if (type == DataType.INT32) return String.valueOf(((int[]) array)[0]);
		if (type == DataType.INT16) return String.valueOf(((short[]) array)[0]);
		if (type == DataType.INT8) return String.valueOf(((byte[]) array)[0]);
		return array.toString();
	}
}