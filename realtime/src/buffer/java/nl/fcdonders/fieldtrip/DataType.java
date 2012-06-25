/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
package nl.fcdonders.fieldtrip;

import java.nio.*;


/** A class for defining FieldTrip data types and routines for
	conversion to Java objects.
*/
public class DataType {
	public static final int UNKNOWN = -1;
	public static final int CHAR    = 0;
	public static final int UINT8   = 1;
	public static final int UINT16  = 2;
	public static final int UINT32  = 3;
	public static final int UINT64  = 4;
	public static final int INT8    = 5;
	public static final int INT16   = 6;
	public static final int INT32   = 7;
	public static final int INT64   = 8;
	public static final int FLOAT32 = 9;
	public static final int FLOAT64 = 10;

	public static final int[] wordSize = {1,1,2,4,8,1,2,4,8,4,8};
	
	public static Object getObject(int type, int numel, ByteBuffer buf) {
		switch(type) {
			case CHAR:
				byte[] strBytes = new byte[numel];
				buf.get(strBytes);
				return new String(strBytes);
			
			case INT8:
			case UINT8:
				byte[] int8array = new byte[numel];
				buf.get(int8array);
				return int8array;
			
			case INT16:
			case UINT16:
				short[] int16array = new short[numel];
				// The following would be faster, but DOES NOT
				// increment the position of the original ByteBuffer!!!
				// buf.asShortBuffer().get(int16array);
				for (int i=0;i<numel;i++) int16array[i] = buf.getShort();
				return int16array;

			case INT32:
			case UINT32:
				int[] int32array = new int[numel];
				for (int i=0;i<numel;i++) int32array[i] = buf.getInt();
				return int32array;
			
			case INT64:
			case UINT64:
				long[] int64array = new long[numel];
				for (int i=0;i<numel;i++) int64array[i] = buf.getLong();
				return int64array;	
			
			case FLOAT32:
				float[] float32array = new float[numel];
				for (int i=0;i<numel;i++) float32array[i] = buf.getFloat();
				return float32array;			
			
			case FLOAT64:
				double[] float64array = new double[numel];
				for (int i=0;i<numel;i++) float64array[i] = buf.getDouble();
				return float64array;			

			default:
				return null;
		}
	}
}	
