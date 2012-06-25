/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
import java.io.*;
import java.nio.*;
import nl.fcdonders.fieldtrip.*;
import java.awt.*;
import java.awt.event.*;

class MarkerGUI {
	Frame frame;
	Button sendButton, connectButton;
	TextField addrField, typeField, valueField;
	BufferClient ftClient;
	Color conColor;

	public MarkerGUI(String address) {
	
		ftClient = new BufferClient();
		
		frame = new Frame("Insert Markers as FieldTrip events");
		frame.setLayout(new GridLayout(4,1));
	
		conColor = new Color(128,255,128);
	
		Panel panA = new Panel();
		Panel panB = new Panel();
		Panel panC = new Panel();
		Panel panD = new Panel();
	
		panA.add(new Label("Address")); 
		addrField = new TextField(address);
		connectButton = new Button("  Connect  ");
		panA.add(addrField);
		panA.add(connectButton);
	
		connectButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				onConnect();
			}
		});
	
	
		panB.add(new Label(" Marker type"));
		typeField = new TextField("Some name", 30); 
		panB.add(typeField);
	
		panC.add(new Label("Marker value"));
		valueField = new TextField("Some value", 30); 
		panC.add(valueField);
	
		sendButton = new Button("Insert marker");
		panD.add(sendButton);
		sendButton.setEnabled(false);
		sendButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				onInsert();
			}
		});
		
		frame.add(panA);
		frame.add(panB);
		frame.add(panC);
		frame.add(panD);
		frame.setSize(400,200);
		frame.setVisible(true);
		
		frame.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
                System.exit(0);
			}
		});
	}
	
	public void disconnect() {
		try {
			ftClient.disconnect();
		}
		catch (IOException e) {}
		connectButton.setLabel("Connect");
		addrField.setBackground(Color.WHITE);
		sendButton.setEnabled(false);
	}
	
	public void connect() {
		try {
			if (ftClient.connect(addrField.getText())) {
				connectButton.setLabel("Disconnect");
				addrField.setBackground(conColor);
				sendButton.setEnabled(true);
			}
		}
		catch (IOException e) {
			disconnect();
		}
	}
	
	public void onConnect() {
		if (ftClient.isConnected()) {
			disconnect();
		} else {
			connect();
		}
	}
	
	public void onInsert() {
		System.out.println("Inserting marker");
		System.out.println(typeField.getText());
		System.out.println(valueField.getText());
		try {
			SamplesEventsCount count = ftClient.poll();
			BufferEvent E = new BufferEvent(typeField.getText(), valueField.getText(), count.nSamples);
			ftClient.putEvent(E);
		}
		catch (IOException e) {
			disconnect();
		}		
	}

	public static void main(String[] args) {
		if (args.length>=1) {
			new MarkerGUI(args[0]);
		} else {
			new MarkerGUI("localhost:1972");
		}
	}
}