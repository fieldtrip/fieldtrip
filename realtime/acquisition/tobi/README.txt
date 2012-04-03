# Overview
The TOBI project defines a standardized interface called TiA (TOBI interface
A) to transmit raw biosignals, supporting multi rate and block-oriented
transmission of different kinds of signals from various acquisition devices
(e.g., EEG, EOG, near-infrared spectroscopy signals, etc.) at the same time
[1].

Interoperability between on streams provided with the TiA and the FieldTrip
buffer is provided with tia2ft. Tia2ft can connect to a TiA server, and serve
the incoming data through a FieldTrip buffer, or optionally, push it to an
external FieldTrip buffer (potentially run on a different device).

Note that only homogeneous signal streams (i.e. with a single sampling rate
and block size) are supported at the moment.


# License
Tia2ft is available under the 3-clause BSD license (see LICENSE.txt). This
license permits commercial use, and is compatible with the GPL license. The
licences of the libraries used in tia2ft are available in their respective
directories.


# Usage


# Compiling

# Linux
To recompile tia2ft, a C++ compiler with the STL and Boost libraries is
required.

To compile tia2ft, you can do the following:

1. build libbuffer.a in /realtime/buffer/src by issuing "make",
2. build tia2ft by issuing "make" in /realtime/acquisition/tobi/.


Then, it can simply be run with:
  $ ./tia2ft.sh --serve-ft-buffer


# Testing tia2ft
To test tia2ft, we can use the TOBI signal server to generate artificial
signals. With the TOBI signal servers version 8ea1376, a TiA serving sine-waves can be started as follows:

  $ ./server.sh bin/server_config.xml

When this sever is running, tia2ft can connect to this server.


# References
[1] IEEE Trans Biomed Eng. 2012 Mar;59(3):852-9. Epub 2011 Dec 21. Proposing a
    standardized protocol for raw biosignal transmission. Breitwieser C, Daly
    I, Neuper C, Müller-Putz GR.

# Related links:
- http://www.ncbi.nlm.nih.gov/pubmed/22194230
- http://dx.doi.org/10.1109/TBME.2011.2174637
- http://www.bcistandards.org/
- http://arxiv.org/pdf/1103.4717.pdf
