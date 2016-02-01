

///////////////////////////////////////////////////////////////////////////////
//
// This class configures and manages the connection to the OpenBCI shield for
// the Arduino.  The connection is implemented via a Serial connection.
// The OpenBCI is configured using single letter text commands sent from the
// PC to the Arduino.  The EEG data streams back from the Arduino to the PC
// continuously (once started).  This class defaults to using binary transfer
// for normal operation.
//
// Created: Chip Audette, Oct 2013
// Modified: through April 2014
// Modified again: Conor Russomanno Sept-Oct 2014
// Modified for Daisy (16-chan) OpenBCI V3: Conor Russomanno Nov 2014
// Modified Daisy Behaviors: Chip Audette Dec 2014
//
// Note: this class now expects the data format produced by OpenBCI V3.
//
/////////////////////////////////////////////////////////////////////////////

//import processing.serial.*;
import java.io.OutputStream; //for logging raw bytes to an output file

final String command_stop = "s";
// final String command_startText = "x";
final String command_startBinary = "b";
final String command_startBinary_wAux = "n";  // already doing this with 'b' now
final String command_startBinary_4chan = "v";  // not necessary now
final String command_activateFilters = "f";  // swithed from 'F' to 'f'  ... but not necessary because taken out of hardware code
final String command_deactivateFilters = "g";  // not necessary anymore 

final String[] command_deactivate_channel = {"1", "2", "3", "4", "5", "6", "7", "8", "q", "w", "e", "r", "t", "y", "u", "i"};
final String[] command_activate_channel = {"!", "@", "#", "$", "%", "^", "&", "*","Q", "W", "E", "R", "T", "Y", "U", "I"};

int channelDeactivateCounter = 0; //used for re-deactivating channels after switching settings...

//everything below is now deprecated...
// final String[] command_activate_leadoffP_channel = {"!", "@", "#", "$", "%", "^", "&", "*"};  //shift + 1-8
// final String[] command_deactivate_leadoffP_channel = {"Q", "W", "E", "R", "T", "Y", "U", "I"};   //letters (plus shift) right below 1-8
// final String[] command_activate_leadoffN_channel = {"A", "S", "D", "F", "G", "H", "J", "K"}; //letters (plus shift) below the letters below 1-8
// final String[] command_deactivate_leadoffN_channel = {"Z", "X", "C", "V", "B", "N", "M", "<"};   //letters (plus shift) below the letters below the letters below 1-8
// final String command_biasAuto = "`";
// final String command_biasFixed = "~";

// ArrayList defaultChannelSettings;

class OpenBCI_ADS1299 {
  
  //final static int DATAMODE_TXT = 0;
  final static int DATAMODE_BIN = 2;
  final static int DATAMODE_BIN_WAUX = 1;  //switched to this value so that receiving Accel data is now the default
  //final static int DATAMODE_BIN_4CHAN = 4;
  
  final static int STATE_NOCOM = 0;
  final static int STATE_COMINIT = 1;
  final static int STATE_SYNCWITHHARDWARE = 2;
  final static int STATE_NORMAL = 3;
  final static int STATE_STOPPED = 4;
  final static int COM_INIT_MSEC = 3000; //you may need to vary this for your computer or your Arduino
  
  //int[] measured_packet_length = {0,0,0,0,0};
  //int measured_packet_length_ind = 0;
  //int known_packet_length_bytes = 0;
  
  final static byte BYTE_START = (byte)0xA0;
  final static byte BYTE_END = (byte)0xC0;
  
  //here is the serial port for this OpenBCI board
  private Serial serial_openBCI = null; 
  private boolean portIsOpen = false;
  
  int prefered_datamode = DATAMODE_BIN_WAUX;
  
  private int state = STATE_NOCOM;
  int dataMode = -1;
  int prevState_millis = 0;
  
  private int nEEGValuesPerPacket = 8; //defined by the data format sent by openBCI boards
  //int nAuxValuesPerPacket = 3; //defined by the data format sent by openBCI boards
  private DataPacket_ADS1299 rawReceivedDataPacket;
  private DataPacket_ADS1299 dataPacket;
  //DataPacket_ADS1299 prevDataPacket;

  private int nAuxValues;
  private boolean isNewDataPacketAvailable = false;
  private OutputStream output; //for debugging  WEA 2014-01-26
  private int prevSampleIndex = 0;
  private int serialErrorCounter = 0;
  
  private final float fs_Hz = 250.0f;  //sample rate used by OpenBCI board...set by its Arduino code
  private final float ADS1299_Vref = 4.5f;  //reference voltage for ADC in ADS1299.  set by its hardware
  private float ADS1299_gain = 24.0;  //assumed gain setting for ADS1299.  set by its Arduino code
  private float scale_fac_uVolts_per_count = ADS1299_Vref / ((float)(pow(2,23)-1)) / ADS1299_gain  * 1000000.f; //ADS1299 datasheet Table 7, confirmed through experiment
  //float LIS3DH_full_scale_G = 4;  // +/- 4G, assumed full scale setting for the accelerometer
  private final float scale_fac_accel_G_per_count = 0.002 / ((float)pow(2,4));  //assume set to +/4G, so 2 mG per digit (datasheet). Account for 4 bits unused
  //final float scale_fac_accel_G_per_count = 1.0;
  private final float leadOffDrive_amps = 6.0e-9;  //6 nA, set by its Arduino code
  
  boolean isBiasAuto = true; //not being used?

  //data related to Conor's setup for V3 boards
  final char[] EOT = {'$','$','$'};
  char[] prev3chars = {'#','#','#'};
  private String defaultChannelSettings = "";
  private String daisyOrNot = ""; 
  private int hardwareSyncStep = 0; //start this at 0...
  private boolean readyToSend = false; //system waits for $$$ after requesting information from OpenBCI board
  private long timeOfLastCommand = 0; //used when sync'ing to hardware

  //some get methods
  public float get_fs_Hz() { return fs_Hz; }
  public float get_Vref() { return ADS1299_Vref; }
  public void set_ADS1299_gain(float _gain) { 
    ADS1299_gain = _gain;  
    scale_fac_uVolts_per_count = ADS1299_Vref / ((float)(pow(2,23)-1)) / ADS1299_gain  * 1000000.0; //ADS1299 datasheet Table 7, confirmed through experiment
  }
  public float get_ADS1299_gain() { return ADS1299_gain; }
  public float get_scale_fac_uVolts_per_count() { return scale_fac_uVolts_per_count; }
  public float get_scale_fac_accel_G_per_count() { return scale_fac_accel_G_per_count; }
  public float get_leadOffDrive_amps() { return leadOffDrive_amps; }
  public String get_defaultChannelSettings() { return defaultChannelSettings; }
  public int get_state() { return state;};
  public boolean get_isNewDataPacketAvailable() { return isNewDataPacketAvailable; }
  
  //constructors
  OpenBCI_ADS1299() {};  //only use this if you simply want access to some of the constants
  OpenBCI_ADS1299(PApplet applet, String comPort, int baud, int nEEGValuesPerOpenBCI, boolean useAux, int nAuxValuesPerPacket) {
    nAuxValues=nAuxValuesPerPacket;
    
    //choose data mode
    println("OpenBCI_ADS1299: prefered_datamode = " + prefered_datamode + ", nValuesPerPacket = " + nEEGValuesPerPacket);
    if (prefered_datamode == DATAMODE_BIN_WAUX) {
      if (!useAux) {
        //must be requesting the aux data, so change the referred data mode
        prefered_datamode = DATAMODE_BIN;
        nAuxValues = 0;
        //println("OpenBCI_ADS1299: nAuxValuesPerPacket = " + nAuxValuesPerPacket + " so setting prefered_datamode to " + prefered_datamode);
      }
    }

    println("OpenBCI_ADS1299: a");

    dataMode = prefered_datamode;

    //allocate space for data packet
    rawReceivedDataPacket = new DataPacket_ADS1299(nEEGValuesPerPacket,nAuxValuesPerPacket);  //this should always be 8 channels
    dataPacket = new DataPacket_ADS1299(nEEGValuesPerOpenBCI,nAuxValuesPerPacket);            //this could be 8 or 16 channels
    //prevDataPacket = new DataPacket_ADS1299(nEEGValuesPerPacket,nAuxValuesPerPacket);
    //set all values to 0 so not null
    for(int i = 0; i < nEEGValuesPerPacket; i++) { 
      rawReceivedDataPacket.values[i] = 0; 
      //prevDataPacket.values[i] = 0; 
    }
    for (int i=0; i < nEEGValuesPerOpenBCI; i++) { dataPacket.values[i]=0; }
    for(int i = 0; i < nAuxValuesPerPacket; i++){
      rawReceivedDataPacket.auxValues[i] = 0;
      dataPacket.auxValues[i] = 0;
      //prevDataPacket.auxValues[i] = 0;
    }

    println("OpenBCI_ADS1299: b");

    //prepare the serial port  ... close if open
    //println("OpenBCI_ADS1299: port is open? ... " + portIsOpen);
    //if(portIsOpen == true) {
    if (isSerialPortOpen()) {
      closeSerialPort();
    }

    println("OpenBCI_ADS1299: i");
    openSerialPort(applet, comPort, baud);
    println("OpenBCI_ADS1299: j");
    
    //open file for raw bytes
    //output = createOutput("rawByteDumpFromProcessing.bin");  //for debugging  WEA 2014-01-26
  }
  
  // //manage the serial port  
  private int openSerialPort(PApplet applet, String comPort, int baud) {
    
    try {
      println("OpenBCI_ADS1299: openSerialPort: attempting to open serial port " + openBCI_portName);
      serial_openBCI = new Serial(applet,comPort,baud); //open the com port
      serial_openBCI.clear(); // clear anything in the com port's buffer    
      portIsOpen = true;
      println("OpenBCI_ADS1299: openSerialPort: port is open (t)? ... " + portIsOpen);
      changeState(STATE_COMINIT);
      return 0;
    } 
    catch (RuntimeException e){
      if (e.getMessage().contains("<init>")) {
        serial_openBCI = null;
        System.out.println("OpenBCI_ADS1299: openSerialPort: port in use, trying again later...");
        portIsOpen = false;
      }
      return 0;
    }
  }

  public int changeState(int newState) {
    state = newState;
    prevState_millis = millis();
    return 0;
  }

  int finalizeCOMINIT() {
    // //wait specified time for COM/serial port to initialize
    // if (state == STATE_COMINIT) {
    //   // println("OpenBCI_ADS1299: finalizeCOMINIT: Initializing Serial: millis() = " + millis());
    //   if ((millis() - prevState_millis) > COM_INIT_MSEC) {
    //     //serial_openBCI.write(command_activates + "\n"); println("Processing: OpenBCI_ADS1299: activating filters");
    //     println("OpenBCI_ADS1299: finalizeCOMINIT: State = NORMAL");
        changeState(STATE_NORMAL);
    //     // startRunning();
    //   }
    // }
    return 0;
  }    
  
  public int closeSDandSerialPort() {
    int returnVal=0;
    
    closeSDFile();
    
    readyToSend = false;
    returnVal = closeSerialPort();
    prevState_millis = 0;  //reset OpenBCI_ADS1299 state clock to use as a conditional for timing at the beginnign of systemUpdate()
    hardwareSyncStep = 0; //reset Hardware Sync step to be ready to go again...
    
    return returnVal;
  }
  
  public int closeSDFile() {
    println("Closing any open SD file. Writing 'j' to OpenBCI.");
    if (serial_openBCI != null) serial_openBCI.write("j"); // tell the SD file to close if one is open...
    delay(100); //make sure 'j' gets sent to the board
    return 0;
  }

  public int closeSerialPort() {
    // if (serial_openBCI != null) {
    println("OpenBCI_ADS1299: closeSerialPort: d");
    portIsOpen = false;
    println("OpenBCI_ADS1299: closeSerialPort: e");
    serial_openBCI.clear();
    println("OpenBCI_ADS1299: closeSerialPort: e2");
    serial_openBCI.stop();
    println("OpenBCI_ADS1299: closeSerialPort: f");
    serial_openBCI = null;
    println("OpenBCI_ADS1299: closeSerialPort: g");
    state = STATE_NOCOM;
    println("OpenBCI_ADS1299: closeSerialPort: h");
    return 0;
  }
  
      
  public void syncWithHardware(int sdSetting){
    switch (hardwareSyncStep) {
      // case 1:
      //   println("openBCI_GUI: syncWithHardware: [0] Sending 'v' to OpenBCI to reset hardware in case of 32bit board...");
      //   serial_openBCI.write('v');
      //   readyToSend = false; //wait for $$$ to iterate... applies to commands expecting a response
      case 1: //send # of channels (8 or 16) ... (regular or daisy setup)
        println("OpenBCI_ADS1299: syncWithHardware: [1] Sending channel count (" + nchan + ") to OpenBCI...");
        if(nchan == 8){
          serial_openBCI.write('c');
        }
        if(nchan == 16){
          serial_openBCI.write('C');
          readyToSend = false;
        }
        break;
      case 2: //reset hardware to default registers 
        println("OpenBCI_ADS1299: syncWithHardware: [2] Reseting OpenBCI registers to default... writing \'d\'...");
        serial_openBCI.write("d"); 
        break;
      case 3: //ask for series of channel setting ASCII values to sync with channel setting interface in GUI
        println("OpenBCI_ADS1299: syncWithHardware: [3] Retrieving OpenBCI's channel settings to sync with GUI... writing \'D\'... waiting for $$$...");
        readyToSend = false; //wait for $$$ to iterate... applies to commands expecting a response
        serial_openBCI.write("D"); 
        break;
      case 4: //check existing registers
        println("OpenBCI_ADS1299: syncWithHardware: [4] Retrieving OpenBCI's full register map for verification... writing \'?\'... waiting for $$$...");
        readyToSend = false; //wait for $$$ to iterate... applies to commands expecting a response
        serial_openBCI.write("?"); 
        break;
      case 5:
        // serial_openBCI.write("j"); // send OpenBCI's 'j' commaned to make sure any already open SD file is closed before opening another one...
        switch (sdSetting){
          case 0: //"Do not write to SD"
            //do nothing
            break;
          case 1: //"5 min max"
            serial_openBCI.write("A");
            break;
          case 2: //"5 min max"
            serial_openBCI.write("S");
            break;
          case 3: //"5 min max"
            serial_openBCI.write("F");
            break;
          case 4: //"5 min max"
            serial_openBCI.write("G");
            break;
          case 5: //"5 min max"
            serial_openBCI.write("H");
            break;
          case 6: //"5 min max"
            serial_openBCI.write("J");
            break;
          case 7: //"5 min max"
            serial_openBCI.write("K");
            break;
          case 8: //"5 min max"
            serial_openBCI.write("L");
            break;
        }
        println("OpenBCI_ADS1299: syncWithHardware: [5] Writing selected SD setting (" + sdSettingString + ") to OpenBCI...");
        if(sdSetting != 0){
          readyToSend = false; //wait for $$$ to iterate... applies to commands expecting a response
        }
        break;
      case 6:
        output("OpenBCI_ADS1299: syncWithHardware: The GUI is done intializing. Click outside of the control panel to interact with the GUI.");
        changeState(STATE_STOPPED);
        systemMode = 10;
        //renitialize GUI if nchan has been updated... needs to be built
        break; 
    }
  }
  
  public void updateSyncState(int sdSetting) {
    //has it been 3000 milliseconds since we initiated the serial port? We want to make sure we wait for the OpenBCI board to finish its setup()
    if ( (millis() - prevState_millis > COM_INIT_MSEC) && (prevState_millis != 0) && (state == openBCI.STATE_COMINIT) ) {
      state = STATE_SYNCWITHHARDWARE;
      timeOfLastCommand = millis();
      serial_openBCI.clear();
      defaultChannelSettings = ""; //clear channel setting string to be reset upon a new Init System
      daisyOrNot = ""; //clear daisyOrNot string to be reset upon a new Init System
      println("OpenBCI_ADS1299: systemUpdate: [0] Sending 'v' to OpenBCI to reset hardware in case of 32bit board...");
      serial_openBCI.write('v');
    }
  
    //if we are in SYNC WITH HARDWARE state ... trigger a command
    if ( (state == STATE_SYNCWITHHARDWARE) && (currentlySyncing == false) ) {
      if(millis() - timeOfLastCommand > 200 && readyToSend == true){
        timeOfLastCommand = millis();
        hardwareSyncStep++;
        syncWithHardware(sdSetting);
      }
    }
  }

  public void sendChar(char val) {
    if (serial_openBCI != null) {
      serial_openBCI.write(key);//send the value as ascii (with a newline character?)
    }
  }
  
  void startDataTransfer(){
    if (serial_openBCI != null) {
      serial_openBCI.clear(); // clear anything in the com port's buffer
      // stopDataTransfer();
      changeState(STATE_NORMAL);  // make sure it's now interpretting as binary
      println("OpenBCI_ADS1299: startDataTransfer(): writing \'" + command_startBinary + "\' to the serial port...");
      serial_openBCI.write(command_startBinary);
    }
  }
  
  public void stopDataTransfer() {
    if (serial_openBCI != null) {
      serial_openBCI.clear(); // clear anything in the com port's buffer
      openBCI.changeState(STATE_STOPPED);  // make sure it's now interpretting as binary
      println("OpenBCI_ADS1299: startDataTransfer(): writing \'" + command_stop + "\' to the serial port...");
      serial_openBCI.write(command_stop);// + "\n");
    }
  }
  
  public boolean isSerialPortOpen() { 
    if (portIsOpen & (serial_openBCI != null)) {
      return true;
    } else {
      return false;
    }
  }
  public boolean isOpenBCISerial(Serial port) {
    if (serial_openBCI == port) {
      return true;
    } else {
      return false;
    }
  }
  
  public void printRegisters() {
    if (serial_openBCI != null) {
      println("OpenBCI_ADS1299: printRegisters(): Writing ? to OpenBCI...");
      openBCI.serial_openBCI.write('?');
    }
  }
  
  //read from the serial port
  public int read() {  return read(false); }
  public int read(boolean echoChar) {
    //println("OpenBCI_ADS1299: read(): State: " + state);
    //get the byte
    byte inByte = byte(serial_openBCI.read());
    
    //write the most recent char to the console
    if (echoChar){  //if not in interpret binary (NORMAL) mode
      // print(".");
      char inASCII = char(inByte); 
      if(isRunning == false && (millis() - timeSinceStopRunning) > 500){
        print(char(inByte));
      }

      //keep track of previous three chars coming from OpenBCI
      prev3chars[0] = prev3chars[1];
      prev3chars[1] = prev3chars[2];
      prev3chars[2] = inASCII;

      if(hardwareSyncStep == 1 && inASCII != '$'){
        daisyOrNot+=inASCII;
        //if hardware returns 8 because daisy is not attached, switch the GUI mode back to 8 channels
        // if(nchan == 16 && char(daisyOrNot.substring(daisyOrNot.length() - 1)) == '8'){
        if(nchan == 16 && daisyOrNot.charAt(daisyOrNot.length() - 1) == '8'){
          verbosePrint(" received from OpenBCI... Switching to nchan = 8 bc daisy is not present...");
          nchan = 8;
        }
      }

      if(hardwareSyncStep == 3 && inASCII != '$'){ //if we're retrieving channel settings from OpenBCI
        defaultChannelSettings+=inASCII;
      }

      //if the last three chars are $$$, it means we are moving on to the next stage of initialization
      if(prev3chars[0] == EOT[0] && prev3chars[1] == EOT[1] && prev3chars[2] == EOT[2]){
        verbosePrint(" > EOT detected...");
        // hardwareSyncStep++;
        prev3chars[2] = '#';
        if(hardwareSyncStep == 3){
          println("OpenBCI_ADS1299: read(): x");
          println(defaultChannelSettings);
          println("OpenBCI_ADS1299: read(): y");
          gui.cc.loadDefaultChannelSettings();
          println("OpenBCI_ADS1299: read(): z");
        }
        readyToSend = true; 
        // println(hardwareSyncStep);
        // syncWithHardware(); //haha, I'm getting very verbose with my naming... it's late...
      }  
    }
    
    //write raw unprocessed bytes to a binary data dump file
    if (output != null) {
      try {
       output.write(inByte);   //for debugging  WEA 2014-01-26
      } catch (IOException e) {
        System.err.println("OpenBCI_ADS1299: read(): Caught IOException: " + e.getMessage());
        //do nothing
      }
    }
    
    interpretBinaryStream(inByte);  //new 2014-02-02 WEA
    return int(inByte);
  }

  /* **** Borrowed from Chris Viegl from his OpenBCI parser for BrainBay
  Modified by Joel Murphy and Conor Russomanno to read OpenBCI data
  Packet Parser for OpenBCI (1-N channel binary format):
  3-byte data values are stored in 'little endian' formant in AVRs
  so this protocol parser expects the lower bytes first.
  Start Indicator: 0xA0
  EXPECTING STANDARD PACKET LENGTH DON'T NEED: Packet_length  : 1 byte  (length = 4 bytes framenumber + 4 bytes per active channel + (optional) 4 bytes for 1 Aux value)
  Framenumber     : 1 byte (Sequential counter of packets)
  Channel 1 data  : 3 bytes 
  ...
  Channel 8 data  : 3 bytes
  Aux Values      : UP TO 6 bytes
  End Indcator    : 0xC0
  TOTAL OF 33 bytes ALL DAY
  ********************************************************************* */
  private int nDataValuesInPacket = 0;
  private int localByteCounter=0;
  private int localChannelCounter=0;
  private int PACKET_readstate = 0;
  // byte[] localByteBuffer = {0,0,0,0};
  private byte[] localAdsByteBuffer = {0,0,0};
  private byte[] localAccelByteBuffer = {0,0};

  void interpretBinaryStream(byte actbyte)  {
    boolean flag_copyRawDataToFullData = false;
    
    //println("OpenBCI_ADS1299: interpretBinaryStream: PACKET_readstate " + PACKET_readstate);
    switch (PACKET_readstate) {
      case 0:  
         //look for header byte  
         if (actbyte == byte(0xA0)) {          // look for start indicator
          // println("OpenBCI_ADS1299: interpretBinaryStream: found 0xA0");
          PACKET_readstate++;
         } 
         break;
      case 1: 
        //check the packet counter
        // println("case 1");
        byte inByte = actbyte;
        rawReceivedDataPacket.sampleIndex = int(inByte); //changed by JAM
        if ((rawReceivedDataPacket.sampleIndex-prevSampleIndex) != 1) {
          if(rawReceivedDataPacket.sampleIndex != 0){  // if we rolled over, don't count as error
            serialErrorCounter++;
            println("OpenBCI_ADS1299: apparent sampleIndex jump from Serial data: " + prevSampleIndex + " to  " + rawReceivedDataPacket.sampleIndex + ".  Keeping packet. (" + serialErrorCounter + ")");
          }
        }
        prevSampleIndex = rawReceivedDataPacket.sampleIndex;
        localByteCounter=0;//prepare for next usage of localByteCounter
        localChannelCounter=0; //prepare for next usage of localChannelCounter
        PACKET_readstate++;
        break;
      case 2: 
        // get ADS channel values 
        // println("case 2");
        localAdsByteBuffer[localByteCounter] = actbyte;
        localByteCounter++;
        if (localByteCounter==3) {
          rawReceivedDataPacket.values[localChannelCounter] = interpret24bitAsInt32(localAdsByteBuffer);
          localChannelCounter++;
          if (localChannelCounter==8) { //nDataValuesInPacket) {  
            // all ADS channels arrived !
            //println("OpenBCI_ADS1299: interpretBinaryStream: localChannelCounter = " + localChannelCounter);
            PACKET_readstate++;
            if (prefered_datamode != DATAMODE_BIN_WAUX) PACKET_readstate++;  //if not using AUX, skip over the next readstate
            localByteCounter = 0;
            localChannelCounter = 0;
            //isNewDataPacketAvailable = true;  //tell the rest of the code that the data packet is complete
          } else { 
            //prepare for next data channel
            localByteCounter=0; //prepare for next usage of localByteCounter
          }
        }
        break;
      case 3:
        // get LIS3DH channel values 2 bytes times 3 axes
        // println("case 3");
        localAccelByteBuffer[localByteCounter] = actbyte;
        localByteCounter++;
        if (localByteCounter==2) {
          rawReceivedDataPacket.auxValues[localChannelCounter]  = interpret16bitAsInt32(localAccelByteBuffer);
          localChannelCounter++;
          if (localChannelCounter==nAuxValues) { //number of accelerometer axis) {  
            // all Accelerometer channels arrived !
            //println("OpenBCI_ADS1299: interpretBinaryStream: Accel Data: " + rawReceivedDataPacket.auxValues[0] + ", " + rawReceivedDataPacket.auxValues[1] + ", " + rawReceivedDataPacket.auxValues[2]);
            PACKET_readstate++;
            localByteCounter = 0;
            //isNewDataPacketAvailable = true;  //tell the rest of the code that the data packet is complete
          } else { 
            //prepare for next data channel
            localByteCounter=0; //prepare for next usage of localByteCounter
          }
        }
        break;
      case 4:
        //look for end byte
        // println("case 4");
        if (actbyte == byte(0xC0)) {    // if correct end delimiter found:
          // println("... 0xC0 found");
          //println("OpenBCI_ADS1299: interpretBinaryStream: found end byte. Setting isNewDataPacketAvailable to TRUE");
          isNewDataPacketAvailable = true; //original place for this.  but why not put it in the previous case block
          flag_copyRawDataToFullData = true;  //time to copy the raw data packet into the full data packet (mainly relevant for 16-chan OpenBCI) 
        } else {
          serialErrorCounter++;
          println("OpenBCI_ADS1299: interpretBinaryStream: Actbyte = " + actbyte);
          println("OpenBCI_ADS1299: interpretBinaryStream: expecteding end-of-packet byte is missing.  Discarding packet. (" + serialErrorCounter + ")");
        }
        PACKET_readstate=0;  // either way, look for next packet
        break;
      default: 
          println("OpenBCI_ADS1299: interpretBinaryStream: Unknown byte: " + actbyte + " .  Continuing...");
          PACKET_readstate=0;  // look for next packet
    }
    
    if (flag_copyRawDataToFullData) {
      copyRawDataToFullData();
    }
  } // end of interpretBinaryStream


  //activate or deactivate an EEG channel...channel counting is zero through nchan-1
  public void changeChannelState(int Ichan,boolean activate) {
    if (serial_openBCI != null) {
      // if ((Ichan >= 0) && (Ichan < command_activate_channel.length)) {
      if ((Ichan >= 0)) {
        if (activate) {
          // serial_openBCI.write(command_activate_channel[Ichan]);
          gui.cc.powerUpChannel(Ichan);
        } else {
          // serial_openBCI.write(command_deactivate_channel[Ichan]);
          gui.cc.powerDownChannel(Ichan);
        }
      }
    }
  }
  
  //deactivate an EEG channel...channel counting is zero through nchan-1
  public void deactivateChannel(int Ichan) {
    if (serial_openBCI != null) {
      if ((Ichan >= 0) && (Ichan < command_deactivate_channel.length)) {
        serial_openBCI.write(command_deactivate_channel[Ichan]);
      }
    }
  }
  
  //activate an EEG channel...channel counting is zero through nchan-1
  public void activateChannel(int Ichan) {
    if (serial_openBCI != null) {
      if ((Ichan >= 0) && (Ichan < command_activate_channel.length)) {
        serial_openBCI.write(command_activate_channel[Ichan]);
      }
    }
  }

  //return the state
  public boolean isStateNormal() { 
    if (state == STATE_NORMAL) { 
      return true;
    } else {
      return false;
    }
  }
  
  // ---- DEPRECATED ---- 
  // public void changeImpedanceState(int Ichan,boolean activate,int code_P_N_Both) {
  //   //println("OpenBCI_ADS1299: changeImpedanceState: Ichan " + Ichan + ", activate " + activate + ", code_P_N_Both " + code_P_N_Both);
  //   if (serial_openBCI != null) {
  //     if ((Ichan >= 0) && (Ichan < command_activate_leadoffP_channel.length)) {
  //       if (activate) {
  //         if ((code_P_N_Both == 0) || (code_P_N_Both == 2)) {
  //           //activate the P channel
  //           serial_openBCI.write(command_activate_leadoffP_channel[Ichan]);
  //         } else if ((code_P_N_Both == 1) || (code_P_N_Both == 2)) {
  //           //activate the N channel
  //           serial_openBCI.write(command_activate_leadoffN_channel[Ichan]);
  //         }
  //       } else {
  //         if ((code_P_N_Both == 0) || (code_P_N_Both == 2)) {
  //           //deactivate the P channel
  //           serial_openBCI.write(command_deactivate_leadoffP_channel[Ichan]);
  //         } else if ((code_P_N_Both == 1) || (code_P_N_Both == 2)) {
  //           //deactivate the N channel
  //           serial_openBCI.write(command_deactivate_leadoffN_channel[Ichan]);
  //         }          
  //       }
  //     }
  //   }
  // }
  
  // public void setBiasAutoState(boolean isAuto) {
  //   if (serial_openBCI != null) {
  //     if (isAuto) {
  //       println("OpenBCI_ADS1299: setBiasAutoState: setting bias to AUTO");
  //       serial_openBCI.write(command_biasAuto);
  //     } else {
  //       println("OpenBCI_ADS1299: setBiasAutoState: setting bias to REF ONLY");
  //       serial_openBCI.write(command_biasFixed);
  //     }
  //   }
  // }
  
  private int interpret24bitAsInt32(byte[] byteArray) {     
    //little endian
    int newInt = ( 
      ((0xFF & byteArray[0]) << 16) |
      ((0xFF & byteArray[1]) << 8) | 
      (0xFF & byteArray[2])
      );
    if ((newInt & 0x00800000) > 0) {
      newInt |= 0xFF000000;
    } else {
      newInt &= 0x00FFFFFF;
    }
    return newInt;
  }
  
  private int interpret16bitAsInt32(byte[] byteArray) {
    int newInt = (
      ((0xFF & byteArray[0]) << 8) |
       (0xFF & byteArray[1])
      );
    if ((newInt & 0x00008000) > 0) {
      newInt |= 0xFFFF0000;
    } else {
      newInt &= 0x0000FFFF;
    }
    return newInt;
  }
  
  
  private int copyRawDataToFullData() {
    //Prior to the 16-chan OpenBCI, we did NOT have rawReceivedDataPacket along with dataPacket...we just had dataPacket.
    //With the 16-chan OpenBCI, where the first 8 channels are sent and then the second 8 channels are sent, we introduced
    //this extra structure so that we could alternate between them.
    //
    //This function here decides how to join the latest data (rawReceivedDataPacket) into the full dataPacket
    
    if (dataPacket.values.length < 2*rawReceivedDataPacket.values.length) {
      //this is an 8 channel board, so simply copy the data
      return rawReceivedDataPacket.copyTo(dataPacket);
    } else {
      //this is 16-channels, so copy the raw data into the correct channels of the new data
      int offsetInd_values = 0;  //this is correct assuming we just recevied a  "board" packet (ie, channels 1-8)
      int offsetInd_aux = 0;     //this is correct assuming we just recevied a  "board" packet (ie, channels 1-8)
      if (rawReceivedDataPacket.sampleIndex % 2 == 0) { // even data packets are from the daisy board
        offsetInd_values = rawReceivedDataPacket.values.length;  //start copying to the 8th slot
        //offsetInd_aux = rawReceivedDataPacket.auxValues.length;  //start copying to the 3rd slot
        offsetInd_aux = 0;  
      }
      return rawReceivedDataPacket.copyTo(dataPacket,offsetInd_values,offsetInd_aux);
    }
  }
  
  public int copyDataPacketTo(DataPacket_ADS1299 target) {
    isNewDataPacketAvailable = false;
    return dataPacket.copyTo(target);
  }
 
 
  private long timeOfLastChannelWrite = 0;
  private int channelWriteCounter = 0;
  private boolean isWritingChannel = false;
  public boolean get_isWritingChannel() { return isWritingChannel; }
  public void configureAllChannelsToDefault() { serial_openBCI.write('d'); };
  public void initChannelWrite(int _numChannel) {  //numChannel counts from zero
      timeOfLastChannelWrite = millis();
      isWritingChannel = true;
  }

  // FULL DISCLAIMER: this method is messy....... very messy... we had to brute force a firmware miscue
  public void writeChannelSettings(int _numChannel,char[][] channelSettingValues) {   //numChannel counts from zero
    if (millis() - timeOfLastChannelWrite >= 50) { //wait 50 milliseconds before sending next character
      verbosePrint("---");
      switch (channelWriteCounter) {
        case 0: //start sequence by send 'x'
          verbosePrint("x" + " :: " + millis());
          serial_openBCI.write('x');
          timeOfLastChannelWrite = millis();
          channelWriteCounter++;
          break;
        case 1: //send channel number
          verbosePrint(str(_numChannel+1) + " :: " + millis());
          if (_numChannel < 8) {
            serial_openBCI.write((char)('0'+(_numChannel+1)));
          }
          if (_numChannel >= 8) {
            //openBCI.serial_openBCI.write((command_activate_channel_daisy[_numChannel-8]));
            serial_openBCI.write((command_activate_channel[_numChannel])); //command_activate_channel holds non-daisy and daisy
          }
          timeOfLastChannelWrite = millis();
          channelWriteCounter++;
          break;
        case 2:  
          verbosePrint(channelSettingValues[_numChannel][channelWriteCounter-2] + " :: " + millis());
          serial_openBCI.write(channelSettingValues[_numChannel][channelWriteCounter-2]);
          //value for ON/OF
          timeOfLastChannelWrite = millis();
          channelWriteCounter++;
          break;
        case 3: 
          verbosePrint(channelSettingValues[_numChannel][channelWriteCounter-2] + " :: " + millis());
          serial_openBCI.write(channelSettingValues[_numChannel][channelWriteCounter-2]);
          //value for ON/OF
          timeOfLastChannelWrite = millis();
          channelWriteCounter++;
          break;
        case 4: 
          verbosePrint(channelSettingValues[_numChannel][channelWriteCounter-2] + " :: " + millis());
          serial_openBCI.write(channelSettingValues[_numChannel][channelWriteCounter-2]);
          //value for ON/OF
          timeOfLastChannelWrite = millis();
          channelWriteCounter++;
          break;
        case 5: 
          verbosePrint(channelSettingValues[_numChannel][channelWriteCounter-2] + " :: " + millis());
          serial_openBCI.write(channelSettingValues[_numChannel][channelWriteCounter-2]);
          //value for ON/OF
          timeOfLastChannelWrite = millis();
          channelWriteCounter++;
          break;
        case 6: 
          verbosePrint(channelSettingValues[_numChannel][channelWriteCounter-2] + " :: " + millis());
          serial_openBCI.write(channelSettingValues[_numChannel][channelWriteCounter-2]);
          //value for ON/OF
          timeOfLastChannelWrite = millis();
          channelWriteCounter++;
          break;
        case 7:
          verbosePrint(channelSettingValues[_numChannel][channelWriteCounter-2] + " :: " + millis());
          serial_openBCI.write(channelSettingValues[_numChannel][channelWriteCounter-2]);
          //value for ON/OF
          timeOfLastChannelWrite = millis();
          channelWriteCounter++;
          break;
        case 8:
          verbosePrint("X" + " :: " + millis());
          serial_openBCI.write('X'); // send 'X' to end message sequence
          timeOfLastChannelWrite = millis();
          channelWriteCounter++;
          break;
        case 9:
          //turn back off channels that were not active before changing channel settings
          switch(channelDeactivateCounter){
            case 0:
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 1:
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 2:
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 3:
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 4:
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 5:
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 6: 
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 7:
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              //check to see if it's 8chan or 16chan ... stop the switch case here if it's 8 chan, otherwise keep going
              if(nchan == 8){
                verbosePrint("done writing channel.");
                isWritingChannel = false;
                channelWriteCounter = 0;
                channelDeactivateCounter = 0;
              } else{
                //keep going
              }
              break;
            case 8:
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 9:
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 10:
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 11:
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 12:
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 13: 
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 14:
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              channelDeactivateCounter++;
              break;
            case 15: 
              if(channelSettingValues[channelDeactivateCounter][0] == '1'){
                verbosePrint("deactivating channel: " + str(channelDeactivateCounter + 1));
                serial_openBCI.write(command_deactivate_channel[channelDeactivateCounter]);
              }
              verbosePrint("done writing channel.");
              isWritingChannel = false;
              channelWriteCounter = 0;
              channelDeactivateCounter = 0;
              break;
          }

          // verbosePrint("done writing channel.");
          // isWritingChannel = false;
          // channelWriteCounter = -1;
          timeOfLastChannelWrite = millis();
          break;
      }
      // timeOfLastChannelWrite = millis();
      // channelWriteCounter++;
    }
  }
  
  private long timeOfLastImpWrite = 0;
  private int impWriteCounter = 0;
  private boolean isWritingImp = false;
  public boolean get_isWritingImp() { return isWritingImp; }
  public void initImpWrite(int _numChannel) {  //numChannel counts from zero
        timeOfLastImpWrite = millis();
        isWritingImp = true;
  }
  public void writeImpedanceSettings(int _numChannel,char[][] impedanceCheckValues) {  //numChannel counts from zero
    //after clicking an impedance button, write the new impedance settings for that channel to OpenBCI
    //after clicking any button, write the new settings for that channel to OpenBCI
    // verbosePrint("Writing impedance settings for channel " + _numChannel + " to OpenBCI!");
    //write setting 1, delay 5ms.. write setting 2, delay 5ms, etc.
    if (millis() - timeOfLastImpWrite >= 50) { //wait 50 milliseconds before sending next character
      verbosePrint("---");
      switch (impWriteCounter) {
        case 0: //start sequence by sending 'z'
          verbosePrint("z" + " :: " + millis());
          serial_openBCI.write('z');
          break;
        case 1: //send channel number
          verbosePrint(str(_numChannel+1) + " :: " + millis());
          if (_numChannel < 8) {
            serial_openBCI.write((char)('0'+(_numChannel+1)));
          }
          if (_numChannel >= 8) {
            //openBCI.serial_openBCI.write((command_activate_channel_daisy[_numChannel-8]));
            serial_openBCI.write((command_activate_channel[_numChannel])); //command_activate_channel holds non-daisy and daisy values
          }
          break;
        case 2: 
        case 3: 
          verbosePrint(impedanceCheckValues[_numChannel][impWriteCounter-2] + " :: " + millis());
          serial_openBCI.write(impedanceCheckValues[_numChannel][impWriteCounter-2]);
          //value for ON/OF
          break;
        case 4:
          verbosePrint("Z" + " :: " + millis());
          serial_openBCI.write('Z'); // send 'X' to end message sequence
          break;
        case 5:
          verbosePrint("done writing imp settings.");
          isWritingImp = false;
          impWriteCounter = -1;
          break;
      }
      timeOfLastImpWrite = millis();
      impWriteCounter++;
    }
  }

 
};  
 
