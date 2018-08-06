package nl.fcdonders.fieldtrip.bufferserver.data;

import nl.fcdonders.fieldtrip.bufferserver.network.NetworkProtocol;
import java.nio.ByteOrder;
import java.nio.ByteBuffer;



/**
 * Wrapper for passing events between DataStore and NetworkProtocol.
 * 
 * @author wieke
 * 
 */
public class Event {
	public final int typeType;
	public final int typeSize;
	public final byte[][] type;
	public final int valueType;
	public final int valueSize;
	public final byte[][] value;
	public final int sample;
	public final int offset;
	public final int duration;
	public final ByteOrder order;

	/**
	 * Constructor
	 * 
	 * @param event
	 *            An event
	 * @param type
	 *            event type in bytes
	 * @param value
	 *            event value in bytes
	 * @param order
	 *            endianess
	 */
	public Event(Event event, byte[][] type, byte[][] value, ByteOrder order) {
		typeType = event.typeType;
		typeSize = event.typeSize;
		valueType = event.valueType;
		valueSize = event.valueSize;
		sample = event.sample;
		offset = event.offset;
		duration = event.duration;
		this.type = type;
		this.value = value;
		this.order = order;
	}

	/**
	 * Constructor
	 * 
	 * @param typeType
	 *            data type of event type
	 * @param typeSize
	 *            number of elements in event type
	 * @param valueType
	 *            data type of event value
	 * @param valueSize
	 *            number of elements in event value
	 * @param sample
	 *            index of related sample
	 * @param offset
	 *            offset in time with regards to sample
	 * @param duration
	 *            duration of the event
	 * @param type
	 *            event type in bytes
	 * @param value
	 *            event value in bytes
	 * @param order
	 *            endianess
	 */
	public Event(int typeType, int typeSize, int valueType, int valueSize,
			int sample, int offset, int duration, byte[][] type,
			byte[][] value, ByteOrder order) {
		this.typeType = typeType;
		this.typeSize = typeSize;
		this.valueType = valueType;
		this.valueSize = valueSize;
		this.sample = sample;
		this.offset = offset;
		this.duration = duration;
		this.type = type;
		this.value = value;
		this.order = order;
	}

	public void serialize(ByteBuffer buf) {
		buf.putInt(typeType);
		buf.putInt(typeSize);
		buf.putInt(valueType);
		buf.putInt(valueSize);
		buf.putInt(sample);
		buf.putInt(offset);
		buf.putInt(duration);
		buf.putInt(typeSize*NetworkProtocol.dataTypeSize(typeType)+valueSize*NetworkProtocol.dataTypeSize(valueType));
		for ( int i=0; i<type.length; i++) buf.put(type[i]);
		for ( int i=0; i<value.length; i++) buf.put(value[i]);
	}	 

}
