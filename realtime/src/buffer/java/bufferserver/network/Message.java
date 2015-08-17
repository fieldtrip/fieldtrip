package nl.fcdonders.fieldtrip.bufferserver.network;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public class Message {
	public final short version;
	public final short type;
	public final ByteBuffer buffer;
	public final ByteOrder order;
	public final long time;

	public Message(final short version, final short type,
			final ByteBuffer buffer, final ByteOrder order, final long time) {
		this.version = version;
		this.type = type;
		this.buffer = buffer;
		this.order = order;
		this.time = time;
	}

	@Override
	public String toString() {
		return "(Version " + Short.toString(version) + ", Type "
				+ Short.toString(type) + ", Size "
				+ Integer.toString(buffer.capacity()) + " time " + time + ")";
	}
}
