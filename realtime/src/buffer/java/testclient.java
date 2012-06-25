/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
import java.io.*;
import java.nio.*;
import nl.fcdonders.fieldtrip.*;

class testclient {
	public static void main(String[] args) throws IOException {
		String hostname = "localhost";
		int port = 1972;
	
		if (args.length>=1) {
			hostname = args[0];
		}
		if (args.length>=2) {
			try {
				port = Integer.parseInt(args[1]);
			}
			catch (NumberFormatException e) {
				port = 0;
			}
			if (port <= 0) {
				System.out.println("Second parameter ("+args[1]+") is not a valid port number.");
				System.exit(1);
			}
		}
		
		BufferClient C = new BufferClient();

		System.out.println("Connecting to "+hostname+":"+port);
		C.connect(hostname, port);
		Header hdr = C.getHeader();
		float[][] data = C.getFloatData(0,hdr.nSamples-1);
		
		
		System.out.println("#channels....: "+hdr.nChans);
		System.out.println("#samples.....: "+hdr.nSamples);
		System.out.println("#events......: "+hdr.nEvents);
		System.out.println("Sampling Freq: "+hdr.fSample);
		System.out.println("data type....: "+hdr.dataType);
		
		for (int n=0;n<hdr.nChans;n++) {
			if (hdr.labels[n] != null) {
				System.out.println("Ch. " + n + ": " + hdr.labels[n]);
			}
		}
				
		
		System.out.println("data:");
		System.out.println("size = " + data.length + " x " + data[0].length);
		//for (int n=0;n<hdr.nSamples;n+=32) {
		//	System.out.println("data[" + n + "][0] = " + data[n][0]);
		//}
		
		/*
		DataDescription descr = new DataDescription();
		ByteBuffer rawBuf = C.getRawData(0,hdr.nSamples-1, descr);
		System.out.println("#samples.....: "+descr.nSamples);
		System.out.println("#channels....: "+descr.nChans);
		System.out.println("data type....: "+descr.dataType);
		System.out.println("sizeBytes....: "+descr.sizeBytes);		
		System.out.println("position.....: "+rawBuf.position());
		System.out.println("capacity.....: "+rawBuf.capacity());
		*/
		
		if (hdr.nEvents > 0) {
			BufferEvent[] evs = C.getEvents(0,hdr.nEvents-1);
			for (int n=0;n<evs.length;n++) {
				System.out.println("Ev: "+n);
				evs[n].print();
			}
		}
		
		/* The following 4 lines have the same effect as the one after this block
		BufferEvent E = new BufferEvent();		
		E.sample = hdr.nSamples-1; // latest sample
		E.setType("Marker"); // string type
		E.setValue(42);      // integer value
		*/
		
		BufferEvent E = new BufferEvent("Marker", 42, hdr.nSamples-1);		
		C.putEvent(E);		
		
		C.disconnect();
	}
}