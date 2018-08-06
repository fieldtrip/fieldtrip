package nl.fcdonders.fieldtrip.bufferserver.data;

public class EventRingBuffer {
	private final Event[] ring;
	private final int capacity;
	private int eventCount = 0;
	private int newPos = 0;

	/**
	 * Constructor
	 * 
	 * @param size
	 *            size of the ring
	 */
	public EventRingBuffer(int size) {
		ring = new Event[size];
		capacity = size;
	}

	/**
	 * Adds an item to the buffer.
	 * 
	 * @param item
	 */
	public synchronized void add(Event item) {
		eventCount++;
		ring[newPos++] = item;

		// If newPos has reached capacity wrap the ring around.
		if (newPos == capacity) {
			newPos = 0;
		}
	}

	/**
	 * Resets the buffer.
	 */
	public synchronized void clear() {
		eventCount = 0;
		newPos = 0;
	}

	/**
	 * Used to get an item from the ring.
	 * 
	 * @param index
	 *            Index ranges from 0 to the number of items added in the ring
	 *            -1.
	 * @return the value at index
	 */
	public synchronized Event get(int index) throws IndexOutOfBoundsException {
		if (index < 0) {
			throw new IndexOutOfBoundsException("Index < 0.");
		}

		if (index < eventCount - capacity) {
			throw new IndexOutOfBoundsException(
					"Index < index of oldest item in buffer.");
		}

		if (index >= eventCount) {
			throw new IndexOutOfBoundsException("Index >= size.");
		}

		if (eventCount < capacity) {
			// Ring hasn't wrapped yet.
			return ring[index];
		} else {
			// Ring has wrapped.

			index -= eventCount - capacity; // Subtract the index of the oldest item
										// still in the ring.

			index += newPos; // Add ring-index of the oldest item still in the
								// ring.

			// Check if index should be wrapped around.
			if (index >= capacity) {
				return ring[index - capacity];
			} else {
				return ring[index];
			}
		}
	}

	/**
	 * Returns the index of the oldest item.
	 * 
	 * @return
	 */
	public synchronized int indexOfOldest() {
		if (eventCount <= capacity) {
			return 0;
		} else {
			return eventCount - capacity;
		}
	}

	/**
	 * Returns the total number of items that have been added to the ring.
	 * 
	 * @return
	 */
	public synchronized int eventCount() {
		return eventCount;
	}

}
