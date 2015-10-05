package nl.fcdonders.fieldtrip.bufferserver.data;

import nl.fcdonders.fieldtrip.bufferserver.exceptions.DataException;
import nl.fcdonders.fieldtrip.bufferserver.network.Request;
import nl.fcdonders.fieldtrip.bufferserver.network.WaitRequest;

public abstract class DataModel {

	/**
	 * Adds a thread, with corresponding request, to the list of listeners of
	 * this dataStore. Once the threshold, as defined in request, had been met
	 * the threads waitOver() function will be called.
	 *
	 * @param thread
	 * @param request
	 */
	public abstract void addWaitRequest(WaitRequest request);

	/**
	 * Removes all data.
	 *
	 * @throws DataException
	 */
	public abstract void flushData() throws DataException;

	/**
	 * Removes all events.
	 *
	 * @throws DataException
	 */
	public abstract void flushEvents() throws DataException;

	/**
	 * Removes the header, and all data & events.
	 *
	 * @throws DataException
	 */
	public abstract void flushHeader() throws DataException;

	/**
	 * Returns all data.
	 *
	 * @return
	 */
	public abstract Data getData() throws DataException;

	/**
	 * Returns the requested data. Throws DataException if impossible.
	 *
	 * @param request
	 *            Start index and end index of the range requested.
	 * @return
	 * @throws DataException
	 */
	public abstract Data getData(Request request) throws DataException;

	/**
	 * Returns the number of events currently stored. Throws DataException if
	 * impossible.
	 *
	 * @return
	 * @throws DataException
	 */
	public abstract int getEventCount() throws DataException;

	/**
	 * Returns all events.
	 *
	 * @return
	 * @throws DataException
	 */
	public abstract Event[] getEvents() throws DataException;

	/**
	 * Returns the requested events. Throws DataException if impossible.
	 *
	 * @param request
	 *            Start index and end index of the range requested.
	 * @return
	 * @throws DataException
	 */
	public abstract Event[] getEvents(Request request) throws DataException;

	/**
	 * Returns the header currently stored. Throws DataException if impossible.
	 *
	 * @return
	 * @throws DataException
	 */
	public abstract Header getHeader() throws DataException;

	/**
	 * Returns the number of currently stored samples. Throws DataException if
	 * impossible.
	 *
	 * @return
	 * @throws DataException
	 */
	public abstract int getSampleCount() throws DataException;

	/**
	 * Returns true if a header exists.
	 *
	 * @return
	 */
	public abstract boolean headerExists();

	/**
	 * Appends the data to the storage. Throws DataException if impossible.
	 *
	 * @param data
	 * @return the number of samples added so far.
	 * @throws DataException
	 */
	public abstract int putData(Data data) throws DataException;

	/**
	 * Appends the events to the storage. Throws DataException if impossible.
	 *
	 * @param events
	 * @return the number of events added so far.
	 * @throws DataException
	 */
	public abstract int putEvents(Event[] events) throws DataException;

	/**
	 * Adds the header to the storage. Throws DataException if impossible.
	 *
	 * @param header
	 * @throws DataException
	 */
	public abstract void putHeader(Header header) throws DataException;


	/**
	 * Cleanup any allocated resources
	 *
	 */
	public abstract void cleanup();

}
