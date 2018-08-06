package nl.fcdonders.fieldtrip.bufferserver.network;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.net.SocketException;
import java.nio.BufferUnderflowException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;

import nl.fcdonders.fieldtrip.bufferserver.data.Chunk;
import nl.fcdonders.fieldtrip.bufferserver.data.Data;
import nl.fcdonders.fieldtrip.bufferserver.data.Event;
import nl.fcdonders.fieldtrip.bufferserver.data.Header;
import nl.fcdonders.fieldtrip.bufferserver.exceptions.ClientException;

/**
 * An impementation of the fieldtrip realtime network protocol. Provides a
 * number of abstract methods that can be used to decode/unwrap incoming
 * communication. And a number of functions that can be used to send outgoing
 * data.
 *
 * @author Wieke Kanters
 *
 */
public class NetworkProtocol {

	/**
	 * Returns the number of bytes in a particular data type.
	 *
	 * @param dataType
	 * @return
	 */
	public static int dataTypeSize(final int dataType) {
		switch (dataType) {
		case CHAR:
		case UINT8:
		case INT8:
			return 1;
		case UINT16:
		case INT16:
			return 2;
		case UINT32:
		case INT32:
		case FLOAT32:
			return 4;
		case UINT64:
		case INT64:
		case FLOAT64:
			return 8;
		}
		return -1;
	}

	/**
	 * Decodes a single extended header chunk from the bytebuffer
	 *
	 * @param buffer
	 * @return
	 * @throws ClientException
	 */
	private static Chunk decodeChunk(final ByteBuffer buffer)
			throws ClientException {
		// Get extended header type
		final int type = buffer.getInt();

		// Get extended header size;
		final int size = buffer.getInt();

		// Check if there are enough bytes remaining
		if (buffer.capacity() - buffer.position() < size) {
			throw new ClientException("Malformed header message.");
		}

		// Grab the remaining bytes in the chunk.

		final byte[] data = new byte[size];
		buffer.get(data);

		return new Chunk(type, size, data);
	}

	/**
	 * Decodes a series of extended header chunks from the bytebuffer.
	 *
	 * @param buffer
	 * @return
	 * @throws ClientException
	 */
	private static Chunk[] decodeChunks(final ByteBuffer buffer)
			throws ClientException {
		final ArrayList<Chunk> chunks = new ArrayList<Chunk>();
		int nChunks = 0;

		// Read events while bytes remain in the buffer.
		try {
			while (buffer.position() < buffer.capacity()) {
				chunks.add(decodeChunk(buffer));
				nChunks++;
			}
		} catch (final BufferUnderflowException e) {
			throw new ClientException("Malformed header message");
		}

		return chunks.toArray(new Chunk[nChunks]);
	}

	/**
	 * Decodes the data from the message. Handles all data as groups of bytes,
	 * does not convert to java primitives.
	 *
	 * @param buf
	 * @return
	 * @throws ClientException
	 */
	public static Data decodeData(final ByteBuffer buffer)
			throws ClientException {
		// Get number of channels
		final int nChans = buffer.getInt();

		// Get number of samples, should be 0
		final int nSamples = buffer.getInt();

		// Get data type
		final int dataType = buffer.getInt();

		// Determine the number of bytes per datapoint.
		final int nBytes = dataTypeSize(dataType);

		// Get size of remaining message.
		final int size = buffer.getInt();

		// Check if size and the number of bytes in the buffer match

		if (buffer.capacity() - buffer.position() != size) {
			throw new ClientException(
					"Defined size of data and actual size do not match.");
		}

		// Check if the number of bytes left in the buffer corresponds to what
		// we expect.
		if (buffer.capacity() - buffer.position() < nSamples * nChans * nBytes) {
			throw new ClientException(
					"Recieved less bytes of data than expected.");
		} else if (buffer.capacity() - buffer.position() > nSamples * nChans
				* nBytes) {
			throw new ClientException(
					"Recieved more bytes of data than expected.");
		}

		// Transfer bytes from the buffer into a nSamples*nChans*nBytes array;
		final byte[][][] data = new byte[nSamples][nChans][nBytes];

		for (int x = 0; x < nSamples; x++) {
			for (int y = 0; y < nChans; y++) {
				for (int z = 0; z < nBytes; z++) {
					data[x][y][z] = buffer.get();
				}
			}
		}

		return new Data(nChans, nSamples, dataType, data, buffer.order());
	}

	/**
	 * Partially decodes a single event from a bytebuffer. Handles type and
	 * value of events as arrays of bytes.
	 *
	 * @param buffer
	 * @return
	 * @throws ClientException
	 */
	private static Event decodeEvent(final ByteBuffer buffer)
			throws ClientException, BufferUnderflowException {
		// Get data type of event type
		final int typeType = buffer.getInt();
		final int typeNBytes = dataTypeSize(typeType);

		// Get number of elements in event type
		final int typeSize = buffer.getInt();

		// Get data type of event value
		final int valueType = buffer.getInt();
		final int valueNBytes = dataTypeSize(valueType);

		// Get number of elements in event value
		final int valueSize = buffer.getInt();

		// Get associated sample
		final int sample = buffer.getInt();

		// Get offset
		final int offset = buffer.getInt();

		// Get duration
		final int duration = buffer.getInt();

		// Get size of remaining data
		final int size = buffer.getInt();

		// Check if size and predicted size are consistent

		if (size != typeSize * typeNBytes + valueSize * valueNBytes) {
			throw new ClientException(
					"Given size and actual size of value and type do not match or malformed event message.");
		}

		if (typeNBytes == -1) {
			throw new ClientException(
					"Wrong type type or malformed event message.");
		}

		if (valueNBytes == -1) {
			throw new ClientException(
					"Wrong value type or malformed event message.");
		}

		// Transfer the remaining bytes in type[][] and value[][]
		final byte[][] type = new byte[typeSize][typeNBytes];

		for (int x = 0; x < typeSize; x++) {
			for (int y = 0; y < typeNBytes; y++) {
				type[x][y] = buffer.get();
			}
		}

		final byte[][] value = new byte[valueSize][valueNBytes];

		for (int x = 0; x < valueSize; x++) {
			for (int y = 0; y < valueNBytes; y++) {
				value[x][y] = buffer.get();
			}
		}

		return new Event(typeType, typeSize, valueType, valueSize, sample,
				offset, duration, type, value, buffer.order());

	}

	/**
	 * Decodes a series of events from the ByteBuffer. Handles event values and
	 * types as bytes.
	 *
	 * @param buffer
	 * @return
	 * @throws ClientException
	 */
	public static Event[] decodeEvents(final ByteBuffer buffer)
			throws ClientException {
		final ArrayList<Event> events = new ArrayList<Event>();
		int nEvents = 0;

		// Read events while bytes remain in the buffer.
		try {
			while (buffer.position() < buffer.capacity()) {
				events.add(decodeEvent(buffer));
				nEvents++;
			}
		} catch (final BufferUnderflowException e) {
			throw new ClientException("Malformed event message");
		}

		return events.toArray(new Event[nEvents]);
	}

	/**
	 * Decodes a header from a bytebuffer.
	 *
	 * @param buf
	 * @return the header object
	 * @throws ClientException
	 *             Thrown if the number of samples/events is higher than 0.
	 */
	public static Header decodeHeader(final ByteBuffer buffer)
			throws ClientException {
		// Get number of channels
		final int nChans = buffer.getInt();

		// Get number of samples, should be 0
		final int nSamples = buffer.getInt();

		if (nSamples != 0) {
			throw new ClientException(
					"Recieved header with more than 0 samples.");
		}

		// Get number of events, should be 0
		final int nEvents = buffer.getInt();

		if (nEvents != 0) {
			throw new ClientException(
					"Recieved header with more than 0 events.");
		}

		// Get sample frequency
		final float fSample = buffer.getFloat();

		// Get data type
		final int dataType = buffer.getInt();

		// Get size of remaining message.
		final int size = buffer.getInt();

		// Check if size matches the remaining bytes

		if (buffer.capacity() - buffer.position() != size) {
			throw new ClientException(
					"Defined size of header chunks and actual size do not match.");
		}

		final Chunk[] chunks = decodeChunks(buffer);

		return new Header(nChans, fSample, dataType, chunks, buffer.order());
	}

	/**
	 * Reads an incoming message and prepares it for further processing.
	 *
	 * @param input
	 * @return A message object containing the version, type and remaining
	 *         bytes. @ * Passed on from input.
	 * @throws ClientException
	 *             Thrown if a version conflict exists between client/server or
	 *             if the client is closing the connection.
	 * @throws IOException
	 */
	public static Message decodeMessage(final BufferedInputStream input)
			throws ClientException, IOException, SocketException {

		// First we determine the endianness of the stream.
		final byte versionByte1 = (byte) input.read();
		final byte versionByte2 = (byte) input.read();

		ByteOrder order;
		if (versionByte1 < versionByte2) {
			order = ByteOrder.BIG_ENDIAN;
		} else {
			order = ByteOrder.LITTLE_ENDIAN;
		}

		// Determine message version
		ByteBuffer buffer = ByteBuffer.allocate(2);
		buffer.order(order);
		buffer.put(versionByte1);
		buffer.put(versionByte2);
		buffer.rewind();
		final short version = buffer.getShort();

		// Check if version corresponds otherwise throw IOException
		if (version == -1) {
			throw new ClientException("Client closing connection.");
		}

		if (version != VERSION) {
			throw new ClientException("Client/Server version conflict, "
					+ "Client Version " + Short.toString(version) + ", "
					+ "Server Version " + Short.toString(VERSION) + ".");
		}

		// Get Message Type
		buffer.rewind();
		loadBuffer(buffer, input, 2);
		final short type = buffer.getShort();

		// Get Message Size
		buffer = ByteBuffer.allocate(4);
		buffer.order(order);
		loadBuffer(buffer, input, 4);
		final int size = buffer.getInt();

		// Get Message body.
		buffer = ByteBuffer.allocate(size);
		buffer.order(order);
		loadBuffer(buffer, input, size);

		return new Message(version, type, buffer, order,
				System.currentTimeMillis());
	}

	/**
	 * Decodes a event/data request.
	 *
	 * @param buf
	 * @return
	 */
	public static Request decodeRequest(final ByteBuffer buffer) {

		// Read begin
		final int begin = buffer.getInt();

		// Read end
		final int end = buffer.getInt();

		return new Request(begin, end);
	}

	/**
	 * Decodes a WaitRequest from the ByteBuffer.
	 *
	 * @param buffer
	 * @return
	 */
	public static WaitRequest decodeWaitRequest(final ByteBuffer buffer) {
		final int nSamples = buffer.getInt();
		final int nEvents = buffer.getInt();
		final int timeout = buffer.getInt();
		return new WaitRequest(nSamples, nEvents, timeout);
	}

	/**
	 * Encodes the Data given the ByteOrder.
	 *
	 * @param output
	 * @param data
	 * @param order
	 *            @
	 */
	public static byte[] encodeData(final Data data, final ByteOrder order) {

		// Create ByteBuffer
		final int nBytes = dataTypeSize(data.dataType);

		final ByteBuffer buffer = ByteBuffer.allocate(8 + 16 + data.size()
				* nBytes);
		buffer.order(order);

		// Add standard message opening
		buffer.putShort(VERSION);
		buffer.putShort(GET_OK);
		buffer.putInt(16 + data.size() * nBytes);

		// Add number of channels
		buffer.putInt(data.nChans);

		// Add number of samples
		buffer.putInt(data.nSamples);

		// Add data type
		buffer.putInt(data.dataType);

		// Add data

		buffer.putInt(data.size() * nBytes);

		final boolean flipOrder = order != data.order && nBytes > 1;

		for (int x = 0; x < data.nSamples; x++) {
			for (int y = 0; y < data.nChans; y++) {
				for (int z = 0; z < nBytes; z++) {
					if (flipOrder) {
						buffer.put(data.data[x][y][nBytes - z - 1]);
					} else {
						buffer.put(data.data[x][y][z]);
					}
				}
			}
		}

		return buffer.array();
	}

	/**
	 * Write an Event to the BufferedOutputStream.
	 *
	 * @param output
	 * @param event
	 * @param order
	 *            @
	 */
	public static byte[] encodeEvents(final Event[] events,
			final ByteOrder order) {

		// Determine total message size
		int totalBufferSize = 8;
		for (final Event event : events) {
			totalBufferSize += 32;
			totalBufferSize += event.typeSize * dataTypeSize(event.typeType);
			totalBufferSize += event.valueSize * dataTypeSize(event.valueType);
		}

		// Create ByteBuffer
		final ByteBuffer buffer = ByteBuffer.allocate(totalBufferSize);
		buffer.order(order);

		// Add standard message opening
		buffer.putShort(VERSION);
		buffer.putShort(GET_OK);
		buffer.putInt(totalBufferSize - 8);

		// Loop through all evens and add them to the buffer.

		for (final Event event : events) {
			// Add event type data type
			buffer.putInt(event.typeType);

			// Add number of elements in event type
			buffer.putInt(event.typeSize);

			// Add event value data type
			buffer.putInt(event.valueType);

			// Add number of elements in event value
			buffer.putInt(event.valueSize);

			// Add associated sample
			buffer.putInt(event.sample);

			// Add offset
			buffer.putInt(event.offset);

			// Add duration
			buffer.putInt(event.duration);

			// Add size of remaining value and type bytes

			final int typeNBytes = dataTypeSize(event.typeType);
			final int valueNBytes = dataTypeSize(event.valueType);

			buffer.putInt(event.typeSize * typeNBytes + event.valueSize
					* valueNBytes);

			// Add type bytes
			boolean flipOrder = order != event.order && typeNBytes > 1;

			for (int x = 0; x < event.typeSize; x++) {
				for (int y = 0; y < typeNBytes; y++) {
					if (flipOrder) {
						buffer.put(event.type[x][typeNBytes - y - 1]);
					} else {
						buffer.put(event.type[x][y]);
					}
				}
			}

			// Add value bytes
			flipOrder = order != event.order && valueNBytes > 1;

			for (int x = 0; x < event.valueSize; x++) {
				for (int y = 0; y < valueNBytes; y++) {
					if (flipOrder) {
						buffer.put(event.value[x][valueNBytes - y - 1]);
					} else {
						buffer.put(event.value[x][y]);
					}
				}
			}
		}

		return buffer.array();
	}

	/**
	 * Encodes a flush error.
	 *
	 * @param order
	 * @return
	 */
	public static byte[] encodeFlushError(final ByteOrder order) {
		final ByteBuffer buffer = ByteBuffer.allocate(8);
		buffer.order(order);

		buffer.putShort(VERSION);
		buffer.putShort(FLUSH_ERR);
		buffer.putInt(0);

		return buffer.array();
	}

	/**
	 * Write a FLUSH_OK to the BufferedOutputStream
	 *
	 * @param output
	 * @param order
	 *            @
	 */
	public static byte[] encodeFlushOkay(final ByteOrder order) {
		final ByteBuffer buffer = ByteBuffer.allocate(8);
		buffer.order(order);

		buffer.putShort(VERSION);
		buffer.putShort(FLUSH_OK);
		buffer.putInt(0);

		return buffer.array();
	}

	/**
	 * Write a GET_ERR to the BufferedOutputStream
	 *
	 * @param output
	 * @param order
	 *            @
	 */
	public static byte[] encodeGetError(final ByteOrder order) {
		final ByteBuffer buffer = ByteBuffer.allocate(8);
		buffer.order(order);

		buffer.putShort(VERSION);
		buffer.putShort(GET_ERR);
		buffer.putInt(0);

		return buffer.array();
	}

	/**
	 * Encodes the Header to the BufferedOutputStream using the given ByteOrder.
	 *
	 * @param output
	 * @param header
	 * @param order
	 *            @
	 */
	public static byte[] encodeHeader(final Header header, final ByteOrder order) {

		// Determine total size
		final int size = 24;

		int chunkSize = 0;

		for (final Chunk chunk : header.chunks) {
			chunkSize += chunk.size + 8;
		}

		// Create a byte buffer.
		final ByteBuffer buffer = ByteBuffer.allocate(8 + size + chunkSize);
		buffer.order(order);

		// Add standard message opening
		buffer.putShort(VERSION);
		buffer.putShort(GET_OK);
		buffer.putInt(24 + chunkSize);

		// Add header information
		buffer.putInt(header.nChans);
		buffer.putInt(header.nSamples);
		buffer.putInt(header.nEvents);
		buffer.putFloat(header.fSample);
		buffer.putInt(header.dataType);
		buffer.putInt(chunkSize);

		// Write extended header chunks

		if (header.nChunks > 0) {
			for (final Chunk chunk : header.chunks) {
				// Write chunk type
				buffer.putInt(chunk.type);

				// Write chunk size
				buffer.putInt(chunk.size);

				// Writ chunk data.
				// In case of Resolutions chunk flip order if necessary.

				final boolean flipOrder = order != header.order;

				if (chunk.type == CHUNK_RESOLUTIONS && flipOrder) {
					for (int i = 0; i < header.nChans; i++) {
						for (int j = 7; j >= 0; j--) {
							buffer.put(chunk.data[i * 8 + j]);
						}
					}
				} else {
					buffer.put(chunk.data);
				}
			}
		}

		return buffer.array();
	}

	/**
	 * Encodes the response to the client for a put error.
	 *
	 * @param output
	 * @param order
	 *            @
	 */
	public static byte[] encodePutError(final ByteOrder order) {
		final ByteBuffer buffer = ByteBuffer.allocate(8);
		buffer.order(order);

		buffer.putShort(VERSION);
		buffer.putShort(PUT_ERR);
		buffer.putInt(0);

		return buffer.array();
	}

	/**
	 * Encodes the response to the client for a successful put.
	 *
	 * @param output
	 * @param order
	 *            @
	 */
	public static byte[] encodePutOkay(final ByteOrder order) {
		final ByteBuffer buffer = ByteBuffer.allocate(8);
		buffer.order(order);

		buffer.putShort(VERSION);
		buffer.putShort(PUT_OK);
		buffer.putInt(0);

		return buffer.array();
	}

	/**
	 * Write a WAIT_ERR to the BufferedOutputStream
	 *
	 * @param output
	 * @param order
	 *            @
	 */
	public static byte[] encodeWaitError(final ByteOrder order) {
		final ByteBuffer buffer = ByteBuffer.allocate(8);
		buffer.order(order);

		buffer.putShort(VERSION);
		buffer.putShort(WAIT_ERR);
		buffer.putInt(0);

		return buffer.array();
	}

	/**
	 * Encodes a WaitRequest given the ByteOrder
	 *
	 * @param output
	 * @param waitRequest
	 * @param order
	 *            @
	 */
	public static byte[] encodeWaitResponse(final int nSamples,
			final int nEvents, final ByteOrder order) {

		// Create ByteBuffer
		final ByteBuffer buffer = ByteBuffer.allocate(16);
		buffer.order(order);

		// Add standard message opening
		buffer.putShort(VERSION);
		buffer.putShort(WAIT_OK);
		buffer.putInt(8);

		// Add nSamples
		buffer.putInt(nSamples);

		// Add nEvents
		buffer.putInt(nEvents);

		return buffer.array();
	}

	/**
	 * Loads a number of bytes from the BufferedInputStream into the ByteBuffer.
	 *
	 * @param buffer
	 * @param input
	 * @param size
	 *            The number of bytes to read. @
	 * @throws IOException
	 */
	private static void loadBuffer(final ByteBuffer buffer,
			final BufferedInputStream input, int size) throws IOException {
		while (size > 0) {
			buffer.put((byte) input.read());
			size--;
		}
		buffer.rewind();
	}

	public static final short VERSION = 1;
	public static final short GET_HDR = 0x201;
	public static final short GET_DAT = 0x202;

	public static final short GET_EVT = 0x203;
	public static final short GET_OK = 0x204;
	public static final short GET_ERR = 0x205;
	public static final short PUT_HDR = 0x101;
	public static final short PUT_DAT = 0x102;

	public static final short PUT_EVT = 0x103;
	public static final short PUT_OK = 0x104;
	public static final short PUT_ERR = 0x105;
	public static final short FLUSH_HDR = 0x301;
	public static final short FLUSH_DAT = 0x302;
	public static final short FLUSH_EVT = 0x303;

	public static final short FLUSH_OK = 0x304;
	public static final short FLUSH_ERR = 0x305;
	public static final short WAIT_DAT = 0x402;
	public static final short WAIT_OK = 0x404;

	public static final short WAIT_ERR = 0x405;
	public static final int CHUNK_UNKNOWN = 0;

	public static final int CHUNK_CHANNEL_NAMES = 1;

	public static final int CHUNK_CHANNEL_FLAGS = 2;

	public static final int CHUNK_RESOLUTIONS = 3;

	public static final int CHUNK_ASCII_KEYVAL = 4;

	public static final int CHUNK_NIFTI1 = 5;

	public static final int CHUNK_SIEMENS_AP = 6;

	public static final int CHUNK_CTF_RES4 = 7;

	public static final int CHUNK_NEUROMAG_HEADER = 8;

	public static final int CHUNK_NEUROMAG_ISOTRAK = 9;

	public static final int CHUNK_NEUROMAG_HPIRESULT = 10;

	public static final int CHAR = 0;

	public static final int UINT8 = 1;

	public static final int UINT16 = 2;

	public static final int UINT32 = 3;

	public static final int UINT64 = 4;

	public static final int INT8 = 5;

	public static final int INT16 = 6;

	public static final int INT32 = 7;

	public static final int INT64 = 8;

	public static final int FLOAT32 = 9;

	public static final int FLOAT64 = 10;

}
