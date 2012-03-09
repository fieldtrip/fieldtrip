
The TCP port numbers for the BrainVision Recorder Remote Data Access (RDA) interface are
 - int16 at ?
 - float at 51244

For debugging a raw TCP dump can be generated with
  $ nc 131.174.xxx.yyy 51244 > dumpfile.bin

To simulate the RDA server, you can run
 $ cat dumpfile.bin | nc -l 51244

More information on netcat (nc) can be found on [1].

[1] http://www.g-loaded.eu/2006/11/06/netcat-a-couple-of-useful-examples/

