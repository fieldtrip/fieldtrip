package nl.fcdonders.fieldtrip.bufferserver;

import java.io.IOException;
import java.net.ServerSocket;
import java.nio.ByteOrder;
import java.util.ArrayList;

import nl.fcdonders.fieldtrip.bufferserver.data.DataModel;
import nl.fcdonders.fieldtrip.bufferserver.data.Header;
import nl.fcdonders.fieldtrip.bufferserver.data.RingDataStore;
import nl.fcdonders.fieldtrip.bufferserver.data.SavingRingDataStore;
import nl.fcdonders.fieldtrip.bufferserver.data.SimpleDataStore;
import nl.fcdonders.fieldtrip.bufferserver.exceptions.DataException;
import nl.fcdonders.fieldtrip.bufferserver.network.ConnectionThread;

/**
 * Buffer class, a thread that opens a serverSocket to listen for connections
 * and starts a connectionThread to handle them.
 *
 * @author wieke
 *
 */
public class BufferServer extends Thread {

	 static final int serverPort =1972;    // default server port
	 static final int dataBufSize=1024*60; // default save samples = 60s @ 1024Hz
	 static final int eventBufSize=10*60;  // default save events  = 60s @ 10Hz

	/**
	 * Main method, starts running a server thread in the current thread.
	 * Handles arguments.
	 *
	 * @param args
	 *            <port> or <port> <nSamplesAndEvents> or <port> <nSamples>
	 *            <nEvents>
	 */
	public static void main(final String[] args) {
		 int logging=0;
		BufferServer buffer=null;
		if (args.length == 0 ){
			 usage();
			 System.exit(1);
		} else if (args.length == 1) {
			 try { 
				  buffer = new BufferServer(Integer.parseInt(args[0]));
			 } catch ( NumberFormatException e ){ // not a number, assume it's the save location
				  buffer = new BufferServer(serverPort, dataBufSize, eventBufSize, args[0]);
			 }
		} else if (args.length == 2) {
			buffer = new BufferServer(Integer.parseInt(args[0]),
					Integer.parseInt(args[1]));
		} else if (args.length == 3) {
			buffer = new BufferServer(Integer.parseInt(args[0]),
					Integer.parseInt(args[1]), Integer.parseInt(args[2]));
		} else if (args.length == 4 ) { // 4rd arg is save location
			 buffer = new BufferServer(Integer.parseInt(args[0]),
												Integer.parseInt(args[1]), Integer.parseInt(args[2]),
												args[3]);
		} else if (args.length == 5 ) { // 5th arg is save location
			 buffer = new BufferServer(Integer.parseInt(args[0]),
												Integer.parseInt(args[1]), Integer.parseInt(args[2]),
												args[3]);
			 logging = Integer.parseInt(args[4]);
		} else {
			 buffer = new BufferServer(serverPort, dataBufSize, eventBufSize);
		}
		buffer.addMonitor(new SystemOutMonitor(logging));
		buffer.run();
		buffer.cleanup();
	}

	 public static void usage(){
		  System.err.println("java -jar BufferServer.jar PORT sampLen eventLen SaveLocation verbosityLevel");
	 }


	private final DataModel dataStore;

	private final int portNumber;
	private ServerSocket serverSocket;
	private boolean disconnectedOnPurpose = false;
	private final ArrayList<ConnectionThread> threads = new ArrayList<ConnectionThread>();
	private FieldtripBufferMonitor monitor = null;
	private int nextClientID = 0;

	/**
	 * Constructor, creates a simple datastore.
	 *
	 * @param portNumber
	 */
	public BufferServer(final int portNumber) {
		this.portNumber = portNumber;
		dataStore = new RingDataStore(dataBufSize,eventBufSize);
		setName("Fieldtrip Buffer Server");
	}

	 public BufferServer(final String path, final int portNumber) {
		this.portNumber = portNumber;
		dataStore = new SavingRingDataStore(dataBufSize,eventBufSize,path);
		setName("Fieldtrip Buffer Server");
	}

	/**
	 * Constructor, creates a ringbuffer that stores nSamplesEvents number of
	 * samples and events.
	 *
	 * @param portNumber
	 * @param nSamplesEvents
	 */
	public BufferServer(final int portNumber, final int nSamplesEvents) {
		this.portNumber = portNumber;
		dataStore = new RingDataStore(nSamplesEvents);
		setName("Fieldtrip Buffer Server");
	}

	/**
	 * Constructor, creates a ringbuffer that stores nSamples of samples and
	 * nEvents of events.
	 *
	 * @param portNumber
	 * @param nSamples
	 * @param nEvents
	 */
	public BufferServer(final int portNumber, final int nSamples, final int nEvents) {
		this.portNumber = portNumber;
		dataStore = new RingDataStore(nSamples, nEvents);
		setName("Fieldtrip Buffer Server");
	}

	public BufferServer(final int portNumber, final int nSamples, final int nEvents, final String path) {
		this.portNumber = portNumber;
		dataStore = new SavingRingDataStore(nSamples, nEvents, path);
		setName("Fieldtrip Buffer Server");
	}

	public void addMonitor(final FieldtripBufferMonitor monitor) {
		this.monitor = monitor;
		for (final ConnectionThread thread : threads) {
			thread.addMonitor(monitor);
		}
	}

	/**
	 * Attempts to close the current serverSocket.
	 *
	 * @throws IOException
	 */
	public void closeConnection() throws IOException {
		serverSocket.close();
	}

	/**
	 * Flushes the events from the datastore.
	 */
	public void flushEvents() {
		try {
			dataStore.flushEvents();
			if (monitor != null) {
				monitor.clientFlushedEvents(-1, System.currentTimeMillis());
			}
		} catch (final DataException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * Flushes the header, data and samples from the dataStore.
	 */
	public void flushHeader() {
		try {
			dataStore.flushHeader();
			if (monitor != null) {
				monitor.clientFlushedHeader(-1, System.currentTimeMillis());
			}
		} catch (final DataException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * Flushes the samples from the datastore.
	 */
	public void flushSamples() {
		try {
			dataStore.flushData();
			if (monitor != null) {
				monitor.clientFlushedData(-1, System.currentTimeMillis());
			}
		} catch (final DataException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * Puts a header into the dataStore.
	 *
	 * @param nChans
	 * @param fSample
	 * @param dataType
	 * @return
	 */
	public boolean putHeader(final int nChans, final float fSample,
			final int dataType) {
		final Header hdr = new Header(nChans, fSample, dataType,
				ByteOrder.nativeOrder());
		try {
			dataStore.putHeader(hdr);
			if (monitor != null) {
				monitor.clientPutHeader(dataType, fSample, nChans, -1,
						System.currentTimeMillis());
			}
			return true;
		} catch (final DataException e) {
			 System.err.println("Error : " + e);
			return false;
		}
	}

	/**
	 * Removes the connection from the list of threads.
	 *
	 * @param connection
	 */
	public synchronized void removeConnection(final ConnectionThread connection) {
		threads.remove(connection);
	}

	/**
	 * Opens a serverSocket and starts listening for connections.
	 */
	@Override
	public void run() {
		 //final Thread mainThread = Thread.currentThread();
		 // Add shutdown hook to force close and flush to disk of open files if interrupted/aborted
		 Runtime.getRuntime().addShutdownHook(new Thread() { public void run() { cleanup(); } });		
		 try {
			serverSocket = new ServerSocket(portNumber);
		} catch (final IOException e) {
			 System.err.println("Could not listen on port " + portNumber);
			 System.err.println(e);
			 e.printStackTrace();
		}
		try {
			 while (true) {
				final ConnectionThread connection = new ConnectionThread(
						nextClientID++, serverSocket.accept(), dataStore, this);
				connection.setName("Fieldtrip Client Thread "
						+ connection.clientAdress);
				connection.addMonitor(monitor);

				synchronized (threads) {
					threads.add(connection);
				}
				connection.start();
			}
		} catch (final IOException e) {
			if (!disconnectedOnPurpose) {
				System.err.println("Server socket disconnected " + portNumber);
				System.err.println(e);
			} else {
				for (final ConnectionThread thread : threads) {
					thread.disconnect();
				}
			}
		}
	}

	/**
	 * Stops the buffer thread and closes all existing client connections.
	 */
	public void stopBuffer() {
		try {
			for (final ConnectionThread thread : threads) {
				thread.disconnect();
			}
			serverSocket.close();
			disconnectedOnPurpose = true;
		} catch (final IOException e) {
		}		
	}

	 public void cleanup(){
		  System.err.println("Running cleanup code");
		  dataStore.cleanup();
	 }
	 
}
