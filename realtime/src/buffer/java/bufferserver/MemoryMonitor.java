package nl.fcdonders.fieldtrip.bufferserver;

/**
 * Simple memory usage monitor. Prints the current memory usage and time to
 * System.out.
 *
 * @author Wieke Kanters
 *
 */
public class MemoryMonitor extends Thread {
	private final Runtime runtime;
	private final BufferServer buffer;

	public MemoryMonitor(final BufferServer buffer) {
		runtime = Runtime.getRuntime();
		this.buffer = buffer;
		setName("Memory Monitor");
	}

	private double Memory() {
		return (runtime.totalMemory() - runtime.freeMemory())
				/ (1024.0 * 1024.0);
	}

	@Override
	public void run() {
		final long start = System.currentTimeMillis();
		double time = start;

		runtime.gc();
		boolean run = true;

		while (run) {
			try {
				sleep(100);
			} catch (final InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			time = (System.currentTimeMillis() - start) / 1000.0;

			System.out.println(time + "\t" + Memory());

			if (time > 10) {
				buffer.stopBuffer();
				run = false;
			}
		}
	}
}