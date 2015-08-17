package nl.fcdonders.fieldtrip.bufferserver.data;
import java.nio.ByteBuffer;

/**
 * Wrapper for the extended header chunks.
 * 
 * @author Wieke Kanters
 * 
 */
public class Chunk {
	public final int type;
	public final int size;
	public final byte[] data;

	/**
	 * Constructor
	 * 
	 * @param type
	 *            The type of header chunk
	 * @param size
	 *            The number of bytes in the chunk
	 * @param data
	 *            The actual data in bytes
	 */
	public Chunk(int type, int size, byte[] data) {
		this.type = type;
		this.size = size;
		this.data = data;
	}

	 public void serialize(ByteBuffer buf){
		  buf.putInt(type);
		  buf.putInt(size);
		  buf.put(data);
	 }
}
