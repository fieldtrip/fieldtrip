# Contents

This dir contains the source files and binaries for audio2ft. Audio2ft grabs
audio data from the sound card, and stream it to a FieldTrip buffer. Further
processing can be performed by clients that connect to the FieldTrip buffer. To
acquire the audio data, the PortAudio library [1] is used.

## TODO:
- Compilation is broken. Perhaps it is designed for an outdated PortAudio
  library.
- Remove config file and replace with getopt.


# Usage

You need to call this tool with a number that selects your sound card and
driver architecture, and the usual command line arguments for selecting the
FieldTrip buffer server address, that is,

  $ audio2ft device [hostname [port]] 
  
where replacing hostname by a minus (-) tells the software to spawn its own
buffer server on the given port. When called without arguments, audio2ft will
list the sound devices found by the PortAudio library, from which you can pick
the desired device for the next call. 

When called with only the device argument, the application will use localhost
with port 1972 as the default settings. Audio data will be captured at a
sampling rate of 44100 Hz, but down sampling can be enabled in the
configuration file config.txt (must be in the same directory), where you can
select channels and attach labels. For more information, please refer to
example_config.txt.


# Compilation 

Audio2ft can be compiled with `make`. Building on Windows is supported through
the MinGW compiler. For windows a pre-compiled DLL is available. 
you might need to compile the buffer library first.

*NOTE THAT COMPILATION SEEMS TO BE BROKEN*. On Ubuntu Linux with PortAudio 1.9
lots of symbols are undefined.


# References

[1] http://www.portaudio.com
