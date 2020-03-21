// This is part of a realtime EEG processing demo. 
//
// The purpose of this device is to simulate an online data stream 
// using a gravity sensor. By shaking the box, the (x, y, z) signals 
// will fluctuate. This is sent over bluetooth using the ModEEG/OpenEEG 
// serial data communication protocol to a computer. On the acquisition 
// computer it can be visualized and processed as if it were normal EEG.
//
// The device is based on a Sparkfun Pro Micro 3V3. 
// The realtime signal is provided by a MMA7361 triple axis accelerometer module.
// The connection to the “EEG acquisition” computer is through a BlueSMiRF bluetooth modem. 
//
// See http://www.fieldtriptoolbox.org/development/realtime/arduino for documentation.

// Copyright (C) 2013, Robert Oostenveld, DCCN
//
// This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
// for the documentation and details.
//
//    FieldTrip is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    FieldTrip is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
//
// $Id$

#define LSB(x) (byte)(x    & 0xff);
#define MSB(x) (byte)(x>>8 & 0xff);

byte hasUSB=0, hasBT=0;
unsigned int x, y, z, tic;
unsigned int xprev=0, yprev=0, zprev=0, count=0;
unsigned long tsample=4; // should be 4 to approximate 256 Hz according to protocol 
byte buf[17];
float alpha = 0.8;

void setup() {

  pinMode(A0, INPUT);  // connected to x 
  pinMode(A1, INPUT);  // connected to y  
  pinMode(A2, INPUT);  // connected to z

  //  Serial.begin(57600);  // connected to USB
  //  tic = millis();
  //  while ((millis()-tic)<1000 && !(Serial));  // max 1 second
  //  hasUSB = (Serial);

  Serial1.begin(57600); // connected to BlueSMiRF
  tic = millis();
  while ((millis()-tic)<1000 && !(Serial1));  // max 1 second
  hasBT = (Serial1);

  buf[ 0] = 0xA5;
  buf[ 1] = 0x5A;
  buf[ 2] = 2; // version
  buf[ 3] = 0; // count
  buf[ 4] = 0; // channel 1
  buf[ 5] = 0;
  buf[ 6] = 0; // channel 2
  buf[ 7] = 0;
  buf[ 8] = 0; // channel 3
  buf[ 9] = 0;
  buf[10] = 0; // channel 4
  buf[11] = 0;
  buf[12] = 0; // channel 5
  buf[13] = 0;
  buf[14] = 0; // channel 6
  buf[15] = 0;
  buf[16] = 0; // switches

  if (hasUSB) {
    Serial.println("[gravity_eeg]");
    Serial.print("hasUSB = ");
    Serial.println(hasUSB);
    Serial.print("hasBT  = ");
    Serial.println(hasBT);
  }
}

void loop() {
  x = analogRead(A0);
  y = analogRead(A1);
  z = analogRead(A2);

  // use an running average with an exponential decay
  x = (1-alpha)*x + alpha*xprev;
  y = (1-alpha)*y + alpha*yprev;
  z = (1-alpha)*z + alpha*zprev;

  // remember the values for the next iteration
  xprev = x;
  yprev = y;
  zprev = z;

  // scale to appropriate values
  x*=32; 
  y*=32; 
  z*=32; 

  // update the data packet
  buf[ 3]++;
  buf[ 4] = MSB(x); // value 1, MSB
  buf[ 5] = LSB(x); // value 1, LSB
  buf[ 6] = MSB(y); // value 2, MSB
  buf[ 7] = LSB(y); // value 2, LSB
  buf[ 8] = MSB(z); // value 3, MSB
  buf[ 9] = LSB(z); // value 3, LSB

  Serial1.flush();
  for (int i=0; i<17; i++) {
    Serial1.write(buf[i]);
  }

  //  Serial1.println(count++);
  //  delay(1000);

  delay(tsample);
}








