package nl.fcdonders.fieldtrip.bufferserver;

public interface FieldtripBufferMonitor {
	public static final int ERROR_PROTOCOL = 0;
	public static final int ERROR_CONNECTION = 1;
	public static final int ERROR_VERSION = 2;

	public void clientClosedConnection(int clientID, long time);

	public void clientContinues(int clientID, long time);

	public void clientError(int clientID, int errorType, long time);

	public void clientFlushedData(int clientID, long time);

	public void clientFlushedEvents(int clientID, long time);

	public void clientFlushedHeader(int clientID, long time);

	public void clientGetEvents(int count, int clientID, long time);

	public void clientGetHeader(int clientID, long time);

	public void clientGetSamples(int count, int clientID, long time);

	public void clientOpenedConnection(int clientID, String adress, long time);

	public void clientPolls(int clientID, long time);

	public void clientPutEvents(int count, int clientID, int diff, long time);

	public void clientPutHeader(int dataType, float fSample, int nChannels,
			int clientID, long time);

	public void clientPutSamples(int count, int clientID, int diff, long time);

	public void clientWaits(int nSamples, int nEvents, int timeout,
			int clientID, long time);
}
