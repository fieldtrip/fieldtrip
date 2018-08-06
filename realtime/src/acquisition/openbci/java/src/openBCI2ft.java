import java.io.*;
import nl.fcdonders.fieldtrip.bufferclient.*;
//import OpenBCI_ADS1299;
//import jssc.SerialPortList;

class openBCI2ft {
	 static int VERB=1; // global verbosity level
	 static int BUFFERSIZE = 65500;

	 //these settings are for a single OpenBCI board
	 static int openBCIbaud = 115200;
	 static int OpenBCI_Nchannels = 8; //normal OpenBCI has 8 channels
	 static int openBCIvaluesperpacket = 8;

	 //use this for when daisy-chaining two OpenBCI boards
	 //int openBCIbaud = 2*115200; //baud rate from the Arduino
	 //final int OpenBCI_Nchannels = 16; //daisy chain has 16 channels

	 //properties of the openBCI board
	 static private float ADS1299_Vref = 4.5f;  //reference voltage for ADC in ADS1299
	 static private float ADS1299_gain = 24;    //assumed gain setting for ADS1299
	 static private float openBCI_impedanceDrive_amps = (float)6.0e-9;  //6 nA
	 boolean isBiasAuto = true;
	 static int n_aux_ifEnabled = 3;  // this is the accelerometer data CHIP 2014-11-03
	 //program constants
	 static char[] defaultSettings={'0','6','0','1','1','0'};


int curDataPacketInd = -1;
int lastReadDataPacketInd = -1;

//related to sync'ing communiction to OpenBCI hardware?
boolean currentlySyncing = false;
long timeOfLastCommand = 0;

	 public static void main(String[] args) throws IOException,InterruptedException {

		  
		  if ( args.length==0 ) {
				System.out.println("openBCI2ft openBCIport bufferhost:bufferport nActiveCh useAux eventDriven buffpacketsize");
		  }
		  
		  // openBCI port
		  String    openBCIport    = null;
		  if (args.length>=1) {
				openBCIport=args[0];
		  }
		  // buffer host:port
		  String buffhostname = "localhost";
		  int buffport = 1972;
		  if (args.length>=2) {
				buffhostname = args[1];
				int sep = buffhostname.indexOf(':');
				if ( sep>0 ) {
					 buffport=Integer.parseInt(buffhostname.substring(sep+1,buffhostname.length()));
					 buffhostname=buffhostname.substring(0,sep);
				}
		  }
		  int nActiveCh=-1;
		  if (args.length>=3) { nActiveCh = Integer.parseInt(args[2]); }		  
		  boolean useAux=true;
		  if (args.length>=4) { 
				useAux = Integer.parseInt(args[3])>0;
				//System.out.println("Warning, samplerate fixed for this hardware. Argument ignored.");
				//sampleRate = Integer.parseInt(args[3]); 
		  }		  
		  int buffdownsample=1;
		  boolean readAll=true; // do we read all data from serial port when some is available
		  boolean serialEvent=false; // do we use event-driven serial comms
		  if (args.length>=5) { 
				serialEvent = Integer.parseInt(args[4])>0;;
		  }
		  // if (args.length>=5) { 
		  // 		//System.out.println("Warning, buff-down-sample currently isn't supported.  Argument ignored.");		 
		  // 		buffdownsample = Integer.parseInt(args[4]); 
				
		  // }		  
		  int buffpacketsize=-1;
		  if (args.length>=6) { buffpacketsize = Integer.parseInt(args[5]); }		  
		  

		  if ( openBCIport == null ) { // list available ports and exit
				System.out.println("No serial port defined.  Current serial ports connected are:");
				String[] portNames = jssc.SerialPortList.getPortNames();
				for(int i = 0; i < portNames.length; i++){
					 System.out.println(portNames[i]);
				}
				System.out.println();
				System.exit(1);
		  }
		  
		  // print the current settings
		  System.out.println("OPENBCI port: " + openBCIport);
		  System.out.println("Buffer server: " + buffhostname + " : " + buffport);

		  // open the connection to the buffer server		  
		  BufferClientClock C = new BufferClientClock();
		  while ( !C.isConnected() ) {
				System.out.println("Connecting to "+buffhostname+":"+buffport);
				try { 
					 C.connect(buffhostname, buffport);
				} catch (IOException ex){
				}
				if ( !C.isConnected() ) { 
					 System.out.println("Couldn't connect. Waiting");
					 Thread.sleep(1000);
				}
		  }
		  System.out.println("Connected to buffer!");
				
		  // open the openBCI port and start the data streaming
		  int nAuxDataValuesPerPacket = n_aux_ifEnabled;
		  int nEEGDataValuesPerPacket = 8;//OpenBCI_ADS1299.get_nChan();
		  OpenBCI_ADS1299 openBCI = null;
		  while ( true ) {
				try {
					 //this also starts the data transfer after XX seconds
					 // BODGE: aux is always sent!!! by board even if we don't want it!
					 openBCI = new OpenBCI_ADS1299(openBCIport,openBCIbaud,nEEGDataValuesPerPacket,true,n_aux_ifEnabled);
					 break;
				} catch (Exception e) {
					 System.out.println("Trying to connect to serial port : " + openBCIport);
					 Thread.sleep(1000); 
				}
		  }
		  System.out.println("Opened the port!");
		  System.out.println("sync with hardware");
		  syncWithHardware(openBCI);
		  System.out.println("Connected to serial port.");

		  // disable unused channels if wanted
		  if ( nActiveCh<0 ) nActiveCh=openBCI.get_nChan();
		  if ( nActiveCh>0 && nActiveCh < openBCI.get_nChan() ) {
				// Channel settings format: is
				// [CHANNEL, POWER_DOWN, GAIN_SET, INPUT_TYPE_SET, BIAS_SET, SRB2_SET, SRB1_SET]
				// **POWER_DOWN**  0 = ON (default),  1 = OFF    
				// **GAIN_SET**    1=2x, 2=4x, 3=6x, 4=8x, 5=12x, 6=24x
				// **INPUT_TYPE_SET**  0=ADSINPUT_NORMAL     	(default),
				//   1=ADSINPUT_SHORTED, 2=ADSINPUT_BIAS_MEAS, 3=ADSINPUT_MVDD, 4=ADSINPUT_TEMP
				//   5=ADSINPUT_TESTSIG, 6=ADSINPUT_BIAS_DRP, 7=ADSINPUT_BIAS_DRN 
				// **BIAS_SET**  0= Remove form BIAS, 1 = Include in BIAS  (default)  
				// **SRB2_SET**  0=Disconnect this input from SRB2, 1=Connect this input to SRB2  (default
				// **SRB1_SET**  0=Disconnect all N inputs from SRB1 (default), 1=Connect all N to SRB1  
				char[][] chSettingsValues=new char[openBCI.get_nChan()][];
				for ( int chi=0; chi<chSettingsValues.length; chi++){
					 //Default: {on, 24x, ADSINPUT_NORMAL, Bias_include, SRB2-include, SRB1-disconnect}
					 chSettingsValues[chi]=new char[defaultSettings.length];
					 for( int ci=0; ci<chSettingsValues[chi].length; ci++) 
						  chSettingsValues[chi][ci]=defaultSettings[ci];
				}
				// disable channels after nActiveCh
				for ( int chi=nActiveCh; chi<openBCI.get_nChan(); chi++){
					 System.out.println("Setting inactive channel : " + chi);
					 chSettingsValues[chi][0]='1'; // power-down
					 //chSettingsValues[chi][2]='5'; // Test-signal
					 openBCI.initChannelWrite(chi);
					 while ( openBCI.get_isWritingChannel() ){ // run the settings loop
						  openBCI.writeChannelSettings(chi,chSettingsValues);
						  Thread.sleep(100); while ( openBCI.read(true,0)>0 );//non-blocking consume output
					 }
				}
				Thread.sleep(100); while ( openBCI.read(true,0)>0 ); // non-blocking consume board output
		  } else {
				System.out.println("Setting channels to default.");
				openBCI.configureAllChannelsToDefault();				
				Thread.sleep(200);
				while ( openBCI.read(true,0)>0 ); // non-blocking consume all board output
		  }
		  System.out.println("sync with hardware");
		  //syncWithHardware(openBCI);
		  System.out.println("Connected to serial port.");

		  // Just in case to ensure we've cleaned the input buffer
		  for ( int i=0; i<10; i++){
				Thread.sleep(100); 
				while ( openBCI.read(true,0)>0 ); // non-blocking consume all board output
		  }
		  // start the data streaming
		  openBCI.startDataTransfer();
		  System.out.println("Data is now running\n-------");
		  
		  // send the header information
		  float sampleRate=openBCI.get_fs_Hz();
		  int   neeg=nActiveCh;
		  int   naux=openBCI.get_nAux();
		  if ( !useAux ) naux=0;
		  int   nch =neeg + naux;
		  Header hdr = new Header(nch,sampleRate,DataType.FLOAT32);
		  if ( VERB>0 ){ System.out.println("Sending header: " + hdr.toString()); }
		  C.putHeader(hdr);

		  if ( !readAll ){
				System.out.println("Disabling readAll");
				openBCI.setReadAllBytes(false);
		  }
		  if ( serialEvent ) {
				System.out.println("Enabling the event listener");
				openBCI.setSerialEventListener(true);
		  } else {
				if ( System.getProperty("os.name").startsWith("Windows") ) {
					 System.out.println("Warning: Blocking mode is unreliable on some Windows versions.");
					 System.out.println("If you are missing samples then switch to event-driven mode.");
				}
		  }

		  if ( buffpacketsize<=0 ) buffpacketsize=(int)Math.ceil(sampleRate/50f);
		  byte[] buffer = new byte[BUFFERSIZE];
		  int openBCIsamp=0;  // current sample number in buffer-packet recieved from openBCI
		  int buffsamp=0; // current sample number in buffer-packet to send to buffer
		  int buffch=0;   // current channel number in buffer-packet to send to buffer
		  int nBlk=0;
		  int nSamp=0;
		  int[] buffpacksz = new int[] {nch,buffpacketsize};
		  float[][] databuff = new float[buffpacketsize][nch];		
		  System.out.println("nCh : " + nActiveCh);
		  System.out.println("#samp/buf : " + buffdownsample + " buff_packet : " + buffpacketsize);
	  		
		  long startT=System.currentTimeMillis(); 
		  long updateT=startT; // time we last printed update
		  long packetT=startT; // time last packet was recieved
		  float packetDur=0; // time between packets in millis
		  // Now do the data forwarding
		  DataPacket_ADS1299 curPacket = 
				new DataPacket_ADS1299(nEEGDataValuesPerPacket,nAuxDataValuesPerPacket); 
		  while ( true ) {			 
				// wait for an OPENBCI message
				// event driven means we don't block... so sleep a sensible amount between
				// samples to not overload the CPU
				while ( !openBCI.get_isNewDataPacketAvailable() ) {
					 if ( VERB>1 ){ 
						  System.out.println((System.currentTimeMillis()-startT)/1000.0 + 
													"Waiting for data packet from openBCI."); 
					 }
					 if ( serialEvent ) { 
						  // BODGE: serial-events and blocking reads don't seem to mix...
						  // so use a non-blocking read to check for new data
						  if ( openBCI.read(false,0) < 0 ){ // non-blocking read new data
								// no data on serial port, so wait for some to come available
								Thread.sleep((int)(.5f*1000f/sampleRate)); // sleep for 50% sample rate
						  }
					 } else {
						  // Use a blocking read to wait for new bytes on the serial port
						  openBCI.read(false,-1);
					 }
				} // got a complete packet				
				nSamp++;
				// Moving average est interpacket duration
				long curT = System.currentTimeMillis();
				packetDur = (System.currentTimeMillis()-packetT)*.01f + packetDur*.99f; 
				packetT   = curT;
				if ( VERB>1 ){ 
					 System.out.println((System.currentTimeMillis()-startT)/1000.0 + 
											  " : Got a data packet from openBCI"); 
				}
				// increment the cursor position
				if ( VERB>0 ){
					 if ( System.currentTimeMillis()-updateT > 2*1000 ) {
						  updateT=System.currentTimeMillis();
						  System.out.println((float)(System.currentTimeMillis()-startT)/1000.0 
													+ "  " + nSamp + "  " + nBlk + "  " + (1000f/packetDur) 
													+ "   (t,samp,blk,Hz)");
					 }
				}
				// get (a copy of) the data just read
				//resets isNewDataPacketAvailable to false
				openBCI.copyDataPacketTo(curPacket); 
				//next, gather the new data into the "little buffer"
				int Ichan=0;
				for (Ichan=0; Ichan < nActiveCh; Ichan++) {   //loop over each cahnnel
					 //scale the data into engineering units ("microvolts") and save to the "little buffer"
					 databuff[buffsamp][Ichan] += curPacket.values[Ichan] * openBCI.get_scale_fac_uVolts_per_count();
				}
				// include the aux data
				for (int auxi=0 ; auxi<naux; auxi++, Ichan++ ){
					 databuff[buffsamp][Ichan] += curPacket.auxValues[auxi];		 
				}
				buffsamp++; // move to next buffer sample

				// move to next buffer sample
				if ( buffsamp >= databuff.length ) { // got a full buffer packet's worth
					 if ( VERB>1 ){ System.out.println("Got buffer packets worth of data. Sending");}
					 // so forward to the buffer
					 // N.B. for efficiency this should probably be double-buffered (sic)
					 C.putData(databuff);
					 // clear out all the old data
					 for ( int i=0; i<databuff.length; i++){
						  for ( int j=0; j<databuff[i].length; j++){
								databuff[i][j]=0;
						  }
					 }
					 buffsamp=0;
					 nBlk++;
				}
		  }
		  // should cleanup correctly... but java doesn't allow unreachable code..
		  // C.disconnect();
	 }
	 
	 static void syncWithHardware(OpenBCI_ADS1299 openBCI) throws InterruptedException{
		  while ( openBCI.get_state() == openBCI.STATE_COMINIT ||
					 openBCI.get_state() == openBCI.STATE_SYNCWITHHARDWARE ){ // start data streaming
				//Argument is SD Card setting: 0= do not write; 1= 5 min; 2= 15 min; 3= 30 min; etc...
				openBCI.updateSyncState(0);
				if ( openBCI.get_state() == openBCI.STATE_COMINIT ) {
					 // Wait for the board to become ready
					 System.out.println("Waiting for device to become ready");
					 Thread.sleep(1000);
				} else if ( !openBCI.get_readyToSend() ) { // read the state update response
					 //System.out.println("**");
					 while ( !openBCI.get_readyToSend() ) {
						  int res= openBCI.read(true,-1);
						  if ( res<0 ) {
								// blocking read the hardware init responses and print them
								System.out.println("Huh?");
						  } else { 
								//System.out.print("[" + res + "]");
						  }
					 }
				} else { // inter-command pause
					 Thread.sleep(100); 
				}
		  }
	 }


} // class
