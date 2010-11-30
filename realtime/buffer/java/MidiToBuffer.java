/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
import java.io.*;
import nl.fcdonders.fieldtrip.*;
import javax.sound.midi.*;

class MidiToBuffer implements Receiver {
	BufferClient ftClient;
	public Transmitter transmitter;

	public MidiToBuffer() {
		ftClient = new BufferClient();
	}
	
	public void disconnect() {
		try {
			ftClient.disconnect();
		}
		catch (IOException e) {}
	}
	
	public boolean connect(String address) {
		try {
			ftClient.connect(address);
		}
		catch (IOException e) {
			System.out.println("Cannot connect to FieldTrip buffer @ " + address);
			return false;
		}
		return true;
	}
	
	public boolean openTransmitter() {
		try {
			transmitter = MidiSystem.getTransmitter();
			transmitter.setReceiver(this);
		}
		catch (MidiUnavailableException e) {
			System.out.println("Exception: " + e);
			transmitter = null;
			return false;
		}
		return true;
	}
	
	public void listDevices() {	
		javax.sound.midi.MidiDevice.Info[] dev = javax.sound.midi.MidiSystem.getMidiDeviceInfo(); 
    	System.out.println("MIDI devices available on this machine:");
    	for (int i = 0; i < dev.length; i++) {
			System.out.print((i+1)+": ");
			System.out.println(dev[i].getName());	
		}
	}
	
	public void close() {
		System.out.println("Closing MIDI receiver.");
	}

	public void send(MidiMessage message, long timeStamp) {
		BufferEvent E = new BufferEvent("MIDI", 0, -1);
		E.setValueUnsigned(message.getMessage());
		System.out.println("\nLength = " + message.getLength());
		System.out.println("Status = " + message.getStatus());
		try {
			// the next three lines can be removed if your buffer 
		    // server runs the auto-sample translation 
			SamplesEventsCount count = ftClient.poll();
			E.sample = count.nSamples;
			System.out.println("Sample = " + E.sample);
		
			ftClient.putEvent(E);
		}
		catch (IOException e) {
			System.out.println("Could not write event.");
		}
	}

	public static void main(String[] args) {
		MidiToBuffer m2b = new MidiToBuffer();
		if (args.length > 0) {
			if (!m2b.connect(args[0])) return;
		}
		m2b.listDevices();
		System.out.println("Trying to open default MIDI IN device...");
		if (!m2b.openTransmitter()) return;
		
		System.out.println("Now listening for MIDI devices. Press q and <enter> to quit.");
		while (true) {
			try {
				int key = System.in.read();
				if (key == 'q') break;
			}
			catch (java.io.IOException e) {}
		}
		System.out.println("Closing...");
		m2b.transmitter.close();
	}
}