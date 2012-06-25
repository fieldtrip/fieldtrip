/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * Simple demo of how to stream audio signals to a Fieldtrip buffer 
 *
 */
import java.io.*;
import nl.fcdonders.fieldtrip.*;
import javax.sound.sampled.*;

class AudioToBuffer {
	BufferClient ftClient;
	TargetDataLine lineIn;

	public AudioToBuffer() {
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
	
	public void listDevices() {
		Mixer.Info[] mixInfo = AudioSystem.getMixerInfo();
    	System.out.println("AUDIO devices available on this machine:");
    	for (int i = 0; i < mixInfo.length; i++) {
			System.out.print((i+1)+": ");
			System.out.println(mixInfo[i]);
			Mixer mixer = AudioSystem.getMixer(mixInfo[i]);
			Line.Info[] lineInfo = mixer.getTargetLineInfo();
			for (int j=0; j < lineInfo.length; j++) {
				System.out.print("   " + (j+1)+": ");
				System.out.println("   " + lineInfo[j]);
			}	
		}
	}
	
	public boolean start() {
		AudioFormat fmt = new AudioFormat(44100.0f, 16, 2, true, false);
		try {
			lineIn = AudioSystem.getTargetDataLine(fmt);
			lineIn.open(fmt);
			lineIn.start();
		}
		catch (LineUnavailableException e) {
			System.out.println(e);
			return false;
		}
		Header hdr = new Header(2, 44100.0f, DataType.INT16);
		try {
			ftClient.putHeader(hdr);
		}
		catch (IOException e) {
			return false;
		}
		return true;
	}
	
	
	public void tick() {
		int na = lineIn.available();
		if (na==0) {
			try {
				Thread.sleep(10);
			}
			catch (InterruptedException e) {}
			return;
		}
		byte[] buf = new byte[na*4];
		lineIn.read(buf, 0, na*4);
		try {
			ftClient.putRawData(na, 2, DataType.INT16, buf);
		}
		catch (IOException e) {
			System.out.println(e);
		}
		System.out.println("Wrote "+na+" samples.");
	}
	
	
	public void stop() {
		lineIn.stop();
	}

	public static void main(String[] args) {
		AudioToBuffer a2b = new AudioToBuffer();
		if (args.length > 0) {
			if (a2b.connect(args[0])==false) return;
		} else {
			System.out.println("Usage:   java AudioToBuffer hostname:port");
			return;
		}
		a2b.listDevices();
		System.out.println("Trying to open default AUDIO IN device...\n");
		if (!a2b.start()) return;
		
		System.out.println("Now streaming audio. Press q and <enter> to quit.\n");
		while (true) {
			a2b.tick();
			try {
				if (System.in.available() > 0) {
					int key = System.in.read();
					if (key == 'q') break;
				}
			}
			catch (java.io.IOException e) {}
		}
		System.out.println("Closing...");
		a2b.stop();
	}
}