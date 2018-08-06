package nl.fcdonders.fieldtrip.bufferserver;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;

public class SystemOutMonitor implements FieldtripBufferMonitor {
	private final HashMap<Integer, String> adresses = new HashMap<Integer, String>();
	private int verbosity=0;
	private int count = 0;
	private final SimpleDateFormat sdf = new SimpleDateFormat(
			"MMM dd,yyyy HH:mm");

	public SystemOutMonitor(int level) {
		adresses.put(-1, "Internal");
	}

	@Override
	public void clientClosedConnection(final int clientID, final long time) {
		 if ( verbosity > 0 )
		System.out.println(sdf.format(new Date(time)) + " Client "
				+ adresses.get(clientID) + " closed connection now " + --count
				+ " connections opened.");
	}

	@Override
	public void clientContinues(final int clientID, final long time) {
		 if ( verbosity > 0 )
		System.out.println(sdf.format(new Date(time)) + " Client "
				+ adresses.get(clientID) + " has continued.");
	}

	@Override
	public void clientError(final int clientID, final int errorType,
			final long time) {
		if (errorType == FieldtripBufferMonitor.ERROR_CONNECTION) {
			System.out.println(sdf.format(new Date(time)) + " Lost client "
					+ adresses.get(clientID) + " connection unexpectidly");
		} else if (errorType == FieldtripBufferMonitor.ERROR_PROTOCOL) {
			System.out.println(sdf.format(new Date(time)) + " Client "
					+ adresses.get(clientID) + " violates protocol");
		} else {
			System.out.println(sdf.format(new Date(time)) + " Client "
					+ adresses.get(clientID) + " has wrong version");
		}
	}

	@Override
	public void clientFlushedData(final int clientID, final long time) {
		 if ( verbosity > 0 )
		System.out.println(sdf.format(new Date(time)) + " Data Flushed by "
				+ adresses.get(clientID));
	}

	@Override
	public void clientFlushedEvents(final int clientID, final long time) {
		 if ( verbosity > 0 )
		System.out.println(sdf.format(new Date(time)) + " Events Flushed by "
				+ adresses.get(clientID));

	}

	@Override
	public void clientFlushedHeader(final int clientID, final long time) {
		 if ( verbosity > 0 )
		System.out.println(sdf.format(new Date(time)) + " Header Flushed by "
				+ adresses.get(clientID));
	}

	@Override
	public void clientGetEvents(final int count, final int clientID,
			final long time) {
		 if ( verbosity > 1 )
		System.out.println(sdf.format(new Date(time)) + " Client "
				+ adresses.get(clientID) + " has been sent " + count
				+ " events.");

	}

	@Override
	public void clientGetHeader(final int clientID, final long time) {
		 if ( verbosity > 1 )
			  System.out.println(sdf.format(new Date(time)) + " Client "
				+ adresses.get(clientID) + " has been sent the header.");
	}

	@Override
	public void clientGetSamples(final int count, final int clientID,
			final long time) {
		 if ( verbosity > 1 )	
			  System.out.println(sdf.format(new Date(time)) + " Client "
				+ adresses.get(clientID) + " has been sent " + count
				+ " samples.");
	}

	@Override
	public void clientOpenedConnection(final int clientID, final String adress,
			final long time) {
		 if ( verbosity > 0 )
			  System.out.println(sdf.format(new Date(time))
				+ " Client opened connection at " + adress + " now " + ++count
				+ " connections opened.");
		adresses.put(clientID, adress);
	}

	@Override
	public void clientPolls(final int clientID, final long time) {
		 if ( verbosity > 1 )	
			  System.out.println(sdf.format(new Date(time)) + " Client "
				+ adresses.get(clientID) + " polls the buffer.");
	}

	@Override
	public void clientPutEvents(final int count, final int clientID,
			final int diff, final long time) {
		 if ( verbosity > 1 )
			  System.out.println(sdf.format(new Date(time)) + " Client "
				+ adresses.get(clientID) + " added " + diff
				+ " events, total now " + count);
	}

	@Override
	public void clientPutHeader(final int dataType, final float fSample,
			final int nChannels, final int clientID, final long time) {
		 if ( verbosity > 0 )
			  System.out.println(sdf.format(new Date(time)) + " Header added by "
				+ adresses.get(clientID) + " datatype " + dataType
				+ " fSample " + fSample + " nChannels " + nChannels);
	}

	@Override
	public void clientPutSamples(final int count, final int clientID,
			final int diff, final long time) {
		 if ( verbosity > 1 )
			  System.out.print(sdf.format(new Date(time)) + " Client "
				+ adresses.get(clientID) + " added " + diff
				+ " samples, total now " + count + "\r");
	}

	@Override
	public void clientWaits(final int nSamples, final int nEvents,
			final int timeout, final int clientID, final long time) {
		 if ( verbosity > 1 )
			  System.out.println(sdf.format(new Date(time)) + " Client "
				+ adresses.get(clientID)
				+ " is now waiting, tresholds are sample count " + nSamples
				+ ", event count " + nEvents + " or timeout " + timeout);
	}

}
