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
	public static final int CHAR    = 0;
	public static final int UINT8   = 1;
	public static final int INT8    = 2;
	public static final int UINT16  = 3;
	public static final int INT16   = 4;
	public static final int UINT32  = 5;
	public static final int INT32   = 6;
	public static final int UINT64  = 7;
	public static final int INT64   = 8;
	public static final int FLOAT32 = 9;
	public static final int FLOAT64 = 10;

	public static final int[] wordSize = {1,1,1,2,2,4,4,8,8,4,8};
	
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
				buf.asShortBuffer().get(int16array);
				return int16array;

			case INT32:
			case UINT32:
				int[] int32array = new int[numel];
				buf.asIntBuffer().get(int32array);
				return int32array;
			
			case INT64:
			case UINT64:
				long[] int64array = new long[numel];
				buf.asLongBuffer().get(int64array);
				return int64array;	

			case FLOAT32:
				float[] float32array = new float[numel];
				buf.asFloatBuffer().get(float32array);
				return float32array;			
			
			case FLOAT64:
				double[] float64array = new double[numel];
				buf.asDoubleBuffer().get(float64array);
				return float64array;			

			default:
				return null;
		}
	}
}