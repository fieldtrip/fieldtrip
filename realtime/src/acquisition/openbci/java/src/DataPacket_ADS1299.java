// Extracted from the dataTypes.pde class.


//////////////////////////////////////
//
// This file contains classes that are helfpul in some way.
//
// Created: Chip Audette, Oct 2013 - Dec 2014
//
/////////////////////////////////////

class DataPacket_ADS1299 {
  int sampleIndex;
  int[] values;
  int[] auxValues;

  //constructor, give it "nValues", which should match the number of values in the
  //data payload in each data packet from the Arduino.  This is likely to be at least
  //the number of EEG channels in the OpenBCI system (ie, 8 channels if a single OpenBCI
  //board) plus whatever auxiliary data the Arduino is sending. 
  DataPacket_ADS1299(int nValues, int nAuxValues) {
    values = new int[nValues];
    auxValues = new int[nAuxValues];
  }
  int printToConsole() {
    System.out.print("printToConsole: DataPacket = ");
    System.out.print(sampleIndex);
    for (int i=0; i < values.length; i++) {
      System.out.print(", " + values[i]);
    }
    for (int i=0; i < auxValues.length; i++) {
      System.out.print(", " + auxValues[i]);
    }
    System.out.println();
    return 0;
  }
  
  int copyTo(DataPacket_ADS1299 target) { return copyTo(target, 0, 0); }
  int copyTo(DataPacket_ADS1299 target, int target_startInd_values, int target_startInd_aux) {
    target.sampleIndex = sampleIndex;
    return copyValuesAndAuxTo(target, target_startInd_values, target_startInd_aux);
  }
  int copyValuesAndAuxTo(DataPacket_ADS1299 target, int target_startInd_values, int target_startInd_aux) {
    int nvalues = values.length;
    for (int i=0; i < nvalues; i++) {
      target.values[target_startInd_values + i] = values[i];
    }
    nvalues = auxValues.length;
    for (int i=0; i < nvalues; i++) {
      target.auxValues[target_startInd_aux + i] = auxValues[i];
    }
    return 0;
  }
};

