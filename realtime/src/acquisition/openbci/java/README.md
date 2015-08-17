This directory contains java based driver software for the openBCI hardware.

This is based extensively on the Processing code from the [openBCI processing GUI](https://github.com/OpenBCI/OpenBCI_Processing), which has been converted from processing-specific-java to general JAVA.

Specificially we have extracted two files from the OpenBCI_GUI processing code:
[OpenBCI_GUI/OpenBCI_ADS1299.pde](https://github.com/OpenBCI/OpenBCI_Processing/tree/master/OpenBCI_GUI) from the directory openBCI_Processing_GUI library.  The orginal file is stored here as OpenBCI_ADS1299.pde.
And [OpenBCI_GUI/dataTypes.pde](https://github.com/OpenBCI/OpenBCI_Processing/tree/master/OpenBCI_GUI) from the directory openBCI_Processing_GUI library.  The orginal file is stored here as dataTypes.pde.

Also to interface between this code and the JSSC serial driver we have extracted the following file from processing itself [src/processing/serial/Serial.java](https://github.com/processing/processing/tree/master/java/libraries/serial/src/processing/serial/).  The orginal source file is stored here as Serial.java.orig

As the openBCI USB dongle presents itself as a serial port this driver uses Java-Simple-Serial-Connection ([JSSC](https://github.com/scream3r/java-simple-serial-connector)) to provide low-level and cross-platform access to the serial ports from java.
