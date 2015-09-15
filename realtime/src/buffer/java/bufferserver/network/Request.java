package nl.fcdonders.fieldtrip.bufferserver.network;

/**
 * Simple wrapper for the details of the get_dat and get_evt messages.
 * 
 * @author Wieke Kanters
 * 
 */
public class Request {
	public final int begin;
	public final int end;

	public Request(int begin, int end) {
		this.begin = begin;
		this.end = end;
	}

}
