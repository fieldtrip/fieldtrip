// This is part of a realtime EEG processing demo. 
//
// The purpose of this particular device is to demonstrate that something can 
// be controlled. It reads the control signal using a wireless RFM12b (433/886 
// MHz) connection and visualizes it with a 10-segment LED array. Instead of 
// driving the LED array, it could also act as a switch or drive a servo motor. 
//
// This device is based on a Arduino Mini R5 at 5V. 
// Wireless connectivity is provided by a RFM12B module.
//
// http://www.fieldtriptoolbox.org/development/realtime/arduino

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

#include <JeeLib.h>

#define LED1 3
#define LED2 4
#define LED3 5
#define LED4 6
#define LED5 7
#define LED6 8
#define LED7 9
#define LED8 A5
#define LED9 A4
// note that A6 and A7 cannot be used for digital out according to http://forum.freetronics.com/viewtopic.php?t=106
#define LED10 A0

#define DURATION 1    // ms, duration for each LED

unsigned int i;
unsigned long latest;
byte prev, val;

void ledarrayon() {
  digitalWrite(LED1, 1);
  digitalWrite(LED2, 1);
  digitalWrite(LED3, 1);
  digitalWrite(LED4, 1);
  digitalWrite(LED5, 1);
  digitalWrite(LED6, 1);
  digitalWrite(LED7, 1);
  digitalWrite(LED8, 1);
  digitalWrite(LED9, 1);
  digitalWrite(LED10, 1);
}

void ledarrayoff() {
  digitalWrite(LED1, 0);
  digitalWrite(LED2, 0);
  digitalWrite(LED3, 0);
  digitalWrite(LED4, 0);
  digitalWrite(LED5, 0);
  digitalWrite(LED6, 0);
  digitalWrite(LED7, 0);
  digitalWrite(LED8, 0);
  digitalWrite(LED9, 0);
  digitalWrite(LED10, 0);
}

/***************************************************************************************************
 * RF12 configuration for 433 Hz 
 **************************************************************************************************/

#define RF12_CS  10
#define RF12_ID  2
#define RF12_GRP 197

/***************************************************************************************************
 * setup 
 **************************************************************************************************/

void setup() {
  Serial.begin(57600);
  while(!Serial) {
    ;
  }
  Serial.println("[gravity_eeg_ledarray]");

  pinMode(LED1, OUTPUT);
  pinMode(LED2, OUTPUT);
  pinMode(LED3, OUTPUT);
  pinMode(LED4, OUTPUT);
  pinMode(LED5, OUTPUT);
  pinMode(LED6, OUTPUT);
  pinMode(LED6, OUTPUT);
  pinMode(LED7, OUTPUT);
  pinMode(LED8, OUTPUT);
  pinMode(LED9, OUTPUT);
  pinMode(LED10, OUTPUT);

  rf12_set_cs(RF12_CS);
  rf12_initialize(RF12_ID, RF12_433MHZ, RF12_GRP);

  Serial.println("RF12 initialized");

  latest = millis();
  val    = 0;
  prev   = 255;
} // setup

int count = 0;
void loop() {

  //  if (Serial.available()) {
  //    val = Serial.read();
  //    latest = millis();
  //  }

  if (rf12_recvDone()) {
    // If the result is true, then a packet has been received and is available for processing. The following global variables will be set:
    // volatile byte rf12_hdr  - Contains the header byte of the received packet - with flag bits and node ID of either the sender or the receiver.
    // volatile byte rf12_len  - The number of data bytes in the packet. A value in the range 0 .. 66.
    // volatile byte rf12_data - A pointer to the received data.
    // volatile byte rf12_crc  - CRC of the received packet, zero indicates correct reception. If != 0 then rf12_hdr, rf12_len, and rf12_data should not be relied upon.
    if (rf12_crc==0 && rf12_len==1) {
      // the packet looks ok
      val    = rf12_data[0];
      latest = millis();

      if (1) {
        //        Serial.println("Received packet");
        //        Serial.print("rf12_len = ");
        //        Serial.println(rf12_len);
        Serial.print("rf12_data = ");
        for (i=0; i<rf12_len; i++){
          Serial.print(rf12_data[i]);
          Serial.print(" ");
        }
      }
      Serial.println();
    }
    else {
      // we only expect one byte, so this seems to be RF junk
      Serial.println("Received junk");
    }      

    if(RF12_WANTS_ACK)
      rf12_sendStart(RF12_ACK_REPLY,0,0);
  }

  if (val!=prev) {
    Serial.print("val = ");
    Serial.println(val);
    prev = val;
  }

  if ((millis()-latest)>2000) {
    // blink every 500 ms
    Serial.println("Timeout");
    float f = millis()/500.0;
    val = ((f - (int)f)<0.5) * 10;
  }

  digitalWrite(LED1, val>0);
  delay(DURATION);
  digitalWrite(LED1, 0);

  digitalWrite(LED2, val>1);
  delay(DURATION);
  digitalWrite(LED2, 0);

  digitalWrite(LED3, val>2);
  delay(DURATION);
  digitalWrite(LED3, 0);

  digitalWrite(LED4, val>3);
  delay(DURATION);
  digitalWrite(LED4, 0);

  digitalWrite(LED5, val>4);
  delay(DURATION);
  digitalWrite(LED5, 0);

  digitalWrite(LED6, val>5);
  delay(DURATION);
  digitalWrite(LED6, 0);

  digitalWrite(LED7, val>6);
  delay(DURATION);
  digitalWrite(LED7, 0);

  digitalWrite(LED8, val>7);
  delay(DURATION);
  digitalWrite(LED8, 0);

  digitalWrite(LED9, val>8);
  delay(DURATION);
  digitalWrite(LED9, 0);

  digitalWrite(LED10, val>9);
  delay(DURATION);
  digitalWrite(LED10, 0);

} // loop









