package nl.fcdonders.fieldtrip.bufferserver.network;

public class WaitRequest {
	public final int nSamples;
	public final int nEvents;
	public final int timeout;

	public WaitRequest(int nSamples, int nEvents, int timeout) {
		 /* Sanity check the inputs... */
		if ( nSamples < 0 ) { this.nSamples=-1; } else { this.nSamples=nSamples; }
		if ( nEvents < 0 )  { this.nEvents=-1;  } else { this.nEvents=nEvents; }
		this.timeout = timeout;
	}

	public synchronized void blockUntilSatisfied(long timeout)
			throws InterruptedException {
		wait(timeout);
	}

	public synchronized void satisfied() {
		notifyAll();
	}

}
