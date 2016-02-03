/* -*- mode: java; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
  PSerial - class for serial port goodness
  Part of the Processing project - http://processing.org

  Copyright (c) 2004-05 Ben Fry & Casey Reas
  Reworked by Gottfried Haider as part of GSOC 2013

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General
  Public License along with this library; if not, write to the
  Free Software Foundation, Inc., 59 Temple Place, Suite 330,
  Boston, MA  02111-1307  USA
*/

//package processing.serial;

//import processing.core.*;

//import java.lang.reflect.*;
import java.util.Map;

import jssc.*;


public class Serial implements SerialPortEventListener {
  //PApplet parent;
  public SerialPort port;
  //Method serialAvailableMethod;
  //Method serialEventMethod;

  byte[] buffer = new byte[32768];
  int inBuffer = 0;
  int readOffset = 0;

  int bufferUntilSize = 1;
  byte bufferUntilByte = 0;

  volatile boolean invokeSerialAvailable = false;
  boolean readAllBytes = true;

  // Things we are currently not exposing:
  // * hardware flow control
  // * state of the RING, RLSD line
  // * sending breaks

  
  public Serial() {
    this("COM1", 9600, 'N', 8, 1);
  }

  
  public Serial(int baudRate) {
    this("COM1", baudRate, 'N', 8, 1);
  }

  
  public Serial(String portName) {
    this(portName, 9600, 'N', 8, 1);
  }

  
  public Serial(String portName, int baudRate) {
    this(portName, baudRate, 'N', 8, 1);
  }

  
  public Serial(String portName, int baudRate, char parity, int dataBits, float stopBits) {
    //this.parent = parent;
    // parent.registerMethod("dispose", this);
    // parent.registerMethod("pre", this);

    // setup parity
    if (parity == 'O') {
      parity = SerialPort.PARITY_ODD;
    } else if (parity == 'E') {
      parity = SerialPort.PARITY_EVEN;
    } else if (parity == 'M') {
      parity = SerialPort.PARITY_MARK;
    } else if (parity == 'S') {
      parity = SerialPort.PARITY_SPACE;
    } else {
      parity = SerialPort.PARITY_NONE;
    }

    // setup stop bits
    int stopBitsIdx = SerialPort.STOPBITS_1;
    if (stopBits == 1.5f) {
      stopBitsIdx = SerialPort.STOPBITS_1_5;
    } else if (stopBits == 2) {
      stopBitsIdx = SerialPort.STOPBITS_2;
    }

    port = new SerialPort(portName);
    try {
      // the native open() call is not using O_NONBLOCK, so this might block for certain operations (see write())
      port.openPort();
      port.setParams(baudRate, dataBits, stopBitsIdx, parity);
      // we could register more events here
    } catch (SerialPortException e) {
      // this used to be a RuntimeException before, so stick with it
      throw new RuntimeException("Error opening serial port " + e.getPortName() + ": " + e.getExceptionType());
    }
  }

  public void addEventListener() throws RuntimeException {
    try {
      // N.B. need to explicitly call this to enable the event-listener interface!
      port.addEventListener(this, SerialPort.MASK_RXCHAR);
      //serialEventMethod = findCallback("serialEvent");
      //serialAvailableMethod = findCallback("serialAvailable");    
    } catch (SerialPortException e) {
      // this used to be a RuntimeException before, so stick with it
      throw new RuntimeException("Error opening serial port " + e.getPortName() + ": " + e.getExceptionType());
    }
  }
  public void removeEventListener() throws RuntimeException {
    try {
      // N.B. need to explicitly call this to enable the event-listener interface!
      port.removeEventListener();
      //serialEventMethod = findCallback("serialEvent");
      //serialAvailableMethod = findCallback("serialAvailable");    
    } catch (SerialPortException e) {
      // this used to be a RuntimeException before, so stick with it
      throw new RuntimeException("Error opening serial port " + e.getPortName() + ": " + e.getExceptionType());
    }
  }

  // private Method findCallback(final String name) {
  //   try {
  //     return parent.getClass().getMethod(name, this.getClass());
  //   } catch (Exception e) {
  //   }
  //   // Permit callback(Object) as alternative to callback(Serial).
  //   try {
  //     return parent.getClass().getMethod(name, Object.class);
  //   } catch (Exception e) {
  //   }
  //   return null;
  // }

  
  public void dispose() {
    stop();
  }

  
  // public void pre() {
  //   if (serialAvailableMethod != null && invokeSerialAvailable) {
  //     invokeSerialAvailable = false;
  //     try {
  //       serialAvailableMethod.invoke(parent, this);
  //     } catch (Exception e) {
  //       System.err.println("Error, disabling serialAvailable() for "+port.getPortName());
  //       System.err.println(e.getLocalizedMessage());
  //       serialAvailableMethod = null;
  //     }
  //   }
  // }


  public int available() {
    return (inBuffer-readOffset);
  }

  
  public void buffer(int size) {
    bufferUntilSize = size;
  }

  
  public void bufferUntil(int inByte) {
    bufferUntilSize = 0;
    bufferUntilByte = (byte)inByte;
  }

  
  public void clear() {
    synchronized (buffer) {
      inBuffer = 0;
      readOffset = 0;
    }
  }


  public boolean getCTS() {
    try {
      return port.isCTS();
    } catch (SerialPortException e) {
      throw new RuntimeException("Error reading the CTS line: " + e.getExceptionType());
    }
  }

  
  public boolean getDSR() {
    try {
      return port.isDSR();
    } catch (SerialPortException e) {
      throw new RuntimeException("Error reading the DSR line: " + e.getExceptionType());
    }
  }

  
  //public static Map<String, String> getProperties(String portName) {
  //  return SerialPortList.getPortProperties(portName);
  //}


  public int last() {
    if (inBuffer == readOffset) {
      return -1;
    }

    synchronized (buffer) {
      int ret = buffer[inBuffer-1] & 0xFF;
      inBuffer = 0;
      readOffset = 0;
      return ret;
    }
  }

  
  public char lastChar() {
    return (char)last();
  }

  
  public static String[] list() {
    // returns list sorted alphabetically, thus cu.* comes before tty.*
    // this was different with RXTX
    return SerialPortList.getPortNames();
  }

  public void setReadAllBytes(boolean val){ synchronized(buffer){readAllBytes=val;} }

	 public int readAll() throws jssc.SerialPortException { 
      // read all available bytes in the serial ports buffer into our local buffer
		  byte[] in=null;
		  synchronized ( port ){ // be sure we don't interleave reads from the port
        synchronized ( buffer ) { // or from the byte-buffer
				// read all bytes in the recieve buffer into our internal buffer
				in = port.readBytes();
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
		  }
		  return -1;
	 }
  
  
	 public int read(int timeout) {
      // read with a (perhaps infinite) time-out
		  int ret=-1;
		  if (inBuffer > readOffset) { // Fast-Path: return data we have in buffer
				synchronized (buffer) {
					 ret = buffer[readOffset++] & 0xFF;
					 if (inBuffer == readOffset) { // reset to start of buffer
						  inBuffer = 0;
						  readOffset = 0;
					 }
				}
		  } else if ( timeout==0 ) { // no-data in buffer, non-blocking check for more
          try { 
             readAll();
          } catch ( jssc.SerialPortException ex ) {
          }
        } else if ( timeout>0 ) { // timeout read
          try{
            synchronized ( port ) {
              byte[] tmp = port.readBytes(1,timeout);						  
              if ( tmp!=null ) ret = tmp[0] & 0xFF;
              if ( readAllBytes ) readAll(); // read the rest if there are any
            }
          } catch ( jssc.SerialPortTimeoutException  ex ) {
            System.err.println("Serial.java: read(timeout) timeout exception");
            ex.printStackTrace();
          } catch ( jssc.SerialPortException ex ) {
            System.err.println("Serial.java: read(timeout) serial port exception");
            ex.printStackTrace();
          }
        } else if ( timeout<0 ) { // block forever read
				try { 
                synchronized ( port ) {
                  byte [] tmp = port.readBytes(1);
                  ret = tmp[0] & 0xFF;
                  if ( readAllBytes ) readAll(); // get the rest if there is any
                }
            } catch ( jssc.SerialPortException ex ) {
					 System.err.println("Serial.java: read(timeout) serial port exception");
					 ex.printStackTrace();
				}
		  }
		  return ret;
	 }

  public int read() {  return read(0);  }

  
  public byte[] readBytes() {
    if (inBuffer == readOffset) {
      return null;
    }

    synchronized (buffer) {
      byte[] ret = new byte[inBuffer-readOffset];
      System.arraycopy(buffer, readOffset, ret, 0, ret.length);
      inBuffer = 0;
      readOffset = 0;
      return ret;
    }
  }

  public int readBytes(byte[] dest) {
    if (inBuffer == readOffset) {
      return 0;
    }

    synchronized (buffer) {
      int toCopy = inBuffer-readOffset;
      if (dest.length < toCopy) {
        toCopy = dest.length;
      }
      System.arraycopy(buffer, readOffset, dest, 0, toCopy);
      readOffset += toCopy;
      if (inBuffer == readOffset) {
        inBuffer = 0;
        readOffset = 0;
      }
      return toCopy;
    }
  }
  

  public byte[] readBytesUntil(int inByte) {
    if (inBuffer == readOffset) {
      return null;
    }

    synchronized (buffer) {
      // look for needle in buffer
      int found = -1;
      for (int i=readOffset; i < inBuffer; i++) {
        if (buffer[i] == (byte)inByte) {
          found = i;
          break;
        }
      }
      if (found == -1) {
        return null;
      }

      int toCopy = found-readOffset+1;
      byte[] dest = new byte[toCopy];
      System.arraycopy(buffer, readOffset, dest, 0, toCopy);
      readOffset += toCopy;
      if (inBuffer == readOffset) {
        inBuffer = 0;
        readOffset = 0;
      }
      return dest;
    }
  }

  
  public int readBytesUntil(int inByte, byte[] dest) {
    if (inBuffer == readOffset) {
      return 0;
    }

    synchronized (buffer) {
      // look for needle in buffer
      int found = -1;
      for (int i=readOffset; i < inBuffer; i++) {
        if (buffer[i] == (byte)inByte) {
          found = i;
          break;
        }
      }
      if (found == -1) {
        return 0;
      }

      // check if bytes to copy fit in dest
      int toCopy = found-readOffset+1;
      if (dest.length < toCopy) {
        System.err.println( "The buffer passed to readBytesUntil() is to small " +
                  "to contain " + toCopy + " bytes up to and including " +
                  "char " + (byte)inByte);
        return -1;
      }
      System.arraycopy(buffer, readOffset, dest, 0, toCopy);
      readOffset += toCopy;
      if (inBuffer == readOffset) {
        inBuffer = 0;
        readOffset = 0;
      }
      return toCopy;
    }
  }

  
  public char readChar() {
    return (char) read();
  }

  
  public String readString() {
    if (inBuffer == readOffset) {
      return null;
    }
    return new String(readBytes());
  }

  
  public String readStringUntil(int inByte) {
    byte temp[] = readBytesUntil(inByte);
    if (temp == null) {
      return null;
    } else {
      return new String(temp);
    }
  }


  public void serialEvent(SerialPortEvent event) {
    if (event.getEventType() == SerialPortEvent.RXCHAR) {
      try {
        readAll();
      } catch (SerialPortException e) {
        throw new RuntimeException("Error reading from serial port " + e.getPortName() + ": " + e.getExceptionType());
      }
    }
  }

  
  public void setDTR(boolean state) {
    // there is no way to influence the behavior of the DTR line when opening the serial port
    // this means that at least on Linux and OS X, Arduino devices are always reset
    try {
      port.setDTR(state);
    } catch (SerialPortException e) {
      throw new RuntimeException("Error setting the DTR line: " + e.getExceptionType());
    }
  }

  
  public void setRTS(boolean state) {
    try {
      port.setRTS(state);
    } catch (SerialPortException e) {
      throw new RuntimeException("Error setting the RTS line: " + e.getExceptionType());
    }
  }

  
  public void stop() {
    try {
      port.closePort();
    } catch (SerialPortException e) {
      // ignored
    }
    inBuffer = 0;
    readOffset = 0;
  }

  
  public void write(byte[] src) {
    try {
      // this might block if the serial device is not yet ready (esp. tty devices under OS X)
      port.writeBytes(src);
      // we used to call flush() here
    } catch (SerialPortException e) {
      throw new RuntimeException("Error writing to serial port " + e.getPortName() + ": " + e.getExceptionType());
    }
  }

  
  public void write(int src) {
    try {
      port.writeInt(src);
    } catch (SerialPortException e) {
      throw new RuntimeException("Error writing to serial port " + e.getPortName() + ": " + e.getExceptionType());
    }
  }

  
  public void write(String src) {
    try {
      port.writeString(src);
    } catch (SerialPortException e) {
      throw new RuntimeException("Error writing to serial port " + e.getPortName() + ": " + e.getExceptionType());
    }
  }
}
