import jssc.SerialPort;
import jssc.SerialPortException;

public class jssc_test implements jssc.SerialPortEventListener {
	 SerialPort port=null;
	 
    public static void main(String[] args) {
		  new jssc_test(args);
    }

	 jssc_test(String[] args){
		  int openBCIPorti=0;
		  String[] portNames = jssc.SerialPortList.getPortNames();
        for(int i = 0; i < portNames.length; i++){
            System.out.println(portNames[i]);
				if ( portNames[i].equals("/dev/ttyUSB0") ) {
					 System.out.println("Found open bci port");
					 openBCIPorti=i;
				}
        }
		  if ( portNames.length==0 ) {
				System.out.println("Couldnt find any serial ports.  Aborting");
				System.exit(-1);
		  }


        port = new SerialPort(portNames[openBCIPorti]);
		  long t0=System.currentTimeMillis();
        try {
            port.openPort();//Open serial port
            port.setParams(jssc.SerialPort.BAUDRATE_115200, 
                                 jssc.SerialPort.DATABITS_8,
                                 jssc.SerialPort.STOPBITS_1,
                                 jssc.SerialPort.PARITY_NONE);
				// N.B. The event listener internally uses a polling loop so takes ~10% CPU
				//port.addEventListener(this,jssc.SerialPort.MASK_RXCHAR);
        }
        catch (SerialPortException ex) {
            System.out.println(ex);
        }

		  System.out.println("openBCI Control Terminal");
		  System.out.println("test commands:");
		  System.out.println("\tb=start streaming,  s=stop streaming");
		  System.out.println("\t?=query settings");
		  System.out.println("\td=set channel defaults");		  
		  System.out.println("\tv=soft-reset the board");		  
		  System.out.println("\tquit=quit this tester program");
		  
		  java.util.Scanner in = new java.util.Scanner(System.in);
        try {
				while (true){
					 // Read a string to send to the port
					 System.out.print("Enter string to send:");
					 String input = in.nextLine();
					 if ( input.equals("quit") ) break;
					 // send to the device
					 port.writeString(input);
					 long tsend = System.currentTimeMillis();
					 // Read the response and print it out
					 //ALT:
					 // port.setParams(115200, 8, 1, 0);//Set params.
					 // Does this block?
					 // random sleep so event callback gets some time to do stuff
					 try{ Thread.sleep(500); } catch ( InterruptedException ex ) { } 
					 while ( System.currentTimeMillis()-tsend < 5000 ) { // 5s between sending
						  long tread=System.currentTimeMillis();
						  // N.B. Blocking read
						  int resp=0; //byte [] resp=null;
						  if ( false ) { // non-blocking read, return null of none to get
								resp = read();
								//resp = port.readBytes(); 
						  } else {
								// Blocking read until there are bytes available
								resp = read(-1);
								//resp = port.readBytes(1); 
						  }
						  //long tret=System.currentTimeMillis();
						  //System.out.print( (tret-t0) + " dt=" + (tret-tread) + " : " );
						  if ( resp>0 ) {
								System.out.print((char)resp);
 							// if ( resp!=null && resp.length>0 ) {
							// 	 for ( int i=0; i<resp.length; i++ ) {
							// 		  System.out.print((char)resp[i]);
							// 	 }
						  } else {
								System.out.print("<null>");
								try{ Thread.sleep(100); } catch ( InterruptedException ex ) { }
						  }
					 }
					 System.out.println();						  
				}
				port.closePort();//Close serial port				
		  } catch ( SerialPortException ex ) {
            System.out.println(ex);				
		  }		  
	 }
	 
	 byte [] buffer=new byte[1024];
	 int inBuffer=0;	 
	 int readOffset=0;
	 public int readAll() throws jssc.SerialPortException { 
		  byte[] in=null;
		  synchronized ( port ){ // be sure we don't interleave reads from the port
				// read all bytes in the recieve buffer into our internal buffer
				in = port.readBytes();
		  }
		  synchronized ( buffer ) {
				if ( in!=null ) {
					 //System.out.println("Got " + in.length + " chars to read.");
					 if ( inBuffer+in.length>buffer.length ) { // increase buffer size
						  // double in size until big enough for all the data
						  int newlen=buffer.length*2; while ( newlen<inBuffer+in.length ) newlen*=2;
						  byte[] temp=new byte[newlen];
						  System.arraycopy(buffer,0,temp,0,inBuffer);
						  buffer=temp;
					 }
					 System.arraycopy(in,0,buffer,inBuffer,in.length);
					 inBuffer+=in.length;
					 return in.length;
				}
		  }
		  return -1;
	 }

	 public int read() { return read(0); }
	 public int read(int timeout) {
		  int ret=-1;
		  if (inBuffer > readOffset) { // return data we have
				synchronized (buffer) {
					 ret = buffer[readOffset++] & 0xFF;
					 if (inBuffer == readOffset) { // reset to start of buffer
						  inBuffer = 0;
						  readOffset = 0;
					 }
				}
		  } else if ( timeout!=0 ) { // wait for some more to read
				try { 
					 byte[] tmp=null;
					 if ( timeout<0 ) { // block for ever read						  
						  System.out.print(".");
						  synchronized ( port ) {
								tmp = port.readBytes(1);
								readAll(); // get the rest if there is any
						  }
					 } else if ( timeout>0 ) { // timeout read
						  //System.out.println("Waiting for bytes to read with timeout");
						  try{
								synchronized ( port ) {
									 tmp = port.readBytes(1,timeout);						  
									 readAll(); // read the rest if there are any
								}
						  } catch ( jssc.SerialPortTimeoutException  ex ) {}
					 }
					 ret = tmp[0] & 0xFF;
				} catch ( jssc.SerialPortException ex ) {
					 System.err.println("Serial.java: read(timeout) serial port exception");
					 ex.printStackTrace();
				}
		  }
		  return ret;
	 }
	 
	 synchronized public void serialEvent(jssc.SerialPortEvent event){
		  System.out.println("Got serial Event" + event.getEventType() + ":" + event.getEventValue());
		  if ( event.getEventType() == jssc.SerialPortEvent.RXCHAR ) {
				synchronized(port){ 
					 synchronized(buffer){
						  try {
								readAll();
						  } catch ( jssc.SerialPortException ex ) {
								System.out.println(ex);
						  }
					 }
				}
		  }
	 }
}
