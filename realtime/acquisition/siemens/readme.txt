Realtime fMRI acquisition tools for Siemens scanners
----------------------------------------------------
(C) 2010, Stefan Klanke
Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

Table of Contents
-----------------
1 - Components and operation principle
2 - Acquiring pixel data on the scanner
3 - Parsing protocol information using "siemensap"
4 - fMRI visualisation clients
5 - Integration in Matlab


1) Components and operation principle
-------------------------------------
The present system for acquiring fMRI data in real-time consists of 
three main blocks: i) a stand-alone executable that runs directly on 
the scanner host, ii) a FieldTrip buffer running on another machine 
in the network, and iii) Matlab scripts or other client applications
that retrieve their data from the buffer. Further to that, a small
modification needs to be made to the applied MR sequence such that
protocol information is written to a specific location on disk.

With the current Siemens scanner software (VB17A), a new file
E:\IMAGE\xx-yyyy\zzzzzzz.PixelData is created on the E:\ drive
of the host computer (the Windows box the scanner is operated from)
immediately after each scan. This file contains pixel data as
unsigned 16-bit integers, where different slices show up as
tiles of a mosaic. The mosaic seems to be always square, and
blank tiles are appended if the number of slices is smaller
than the number of tiles in the mosaic. 

Example: The MR sequence is set up to scan N=32 slices with
readout resolution R=64 pixels and phase resolution P=48 pixels.
In this case, the mosaic will contain 6x6 tiles with 4 empty
tiles marked by "--" and slices ordered as follows:

01  02  03  04  05  06
07  08  09  10  11  12
13  14  15  16  17  18
19  20  21  22  23  24
25  26  27  28  29  30
31  32  --  --  --  --

The pixel dimensions of the mosaic will be (64*6) x (48*6),
that is, 384 x 288, and thus the total number of pixels is
110592, corresponding to a file size of 221184 bytes. Within
the file, the pixels are written row after row, that is,
the first 768 bytes contain the 384 pixels of the first row,
corresponding to the first rows of slices 01-06, and so on.

Within the FieldTrip buffer, each scan is represented as ONE
sample with RxPxN channels, with data ordering as in Matlab, 
that is, the pixel data is reshaped such that slices (and their 
rows) are contiguous in memory, and empty tiles are dropped. 
For the above example, we would have 64x48x32 = 98304 "channels".
The data format is kept as INT16_T.

Along with each sample, the acquisition tool also writes a "timestamp" 
event to the buffer, which is represented by seconds.microseconds 
since 1970 (UNIX time). Note that this timestamp refers to the
system time (of the host) at the moment the PixelData is read from
disk, which is a few milliseconds later than the actual scan time.
The true latency involved here has not yet been accurately measured
(how?). The software components involved here are described in
section 2.

At the start of the MR sequence, protocol information is written 
into the header of the FieldTrip buffer in the format of Siemens
ASCII protocol (.pro) files. A library for parsing these files
is provided (section 3).

Once the data is in the FieldTrip buffer, all the usual methods
for access can be used, e.g., the read_data Matlab function. Two
demos for accessing the data from C++ are described in section 4.


2) Acquiring pixel data on the scanner
--------------------------------------
In order to react efficiently to a new scan, which shows up as
a .PixelData file with name and location not known in advance, the 
Windows API function ReadDirectoryChangesW is used to monitor E:\IMAGE 
and all of its subdirectories. Whenever a new file is created or 
modified anywhere in that tree, a Windows event is triggered and 
the corresponding path is made available. This mechanism is 
wrapped up in the C++ class FolderWatcher.

A second C++ class, PixelDataGrabber, encapsulates the actual real-time
fMRI acquisition mechanism based on the FolderWatcher and client-side 
code of the FieldTrip buffer. Detailed Doxygen-style documentation is 
provided in PixelDataGrabber.h, and developers can also look at 
pixeldata_to_remote_buffer.cc for a simple example of using this
class in a command-line application. A slightly more complex program
is compiled from gui_streamer.cc, which combines the PixelDataGrabber
with a small GUI written with FLTK (http://www.fltk.org). This program
provides a few buttons for starting and stopping to monitor for new
files, and to connect to/disconnect from a FieldTrip buffer.

How does the PixelDataGrabber determine the number of slices and
their dimensions? For this to work, the best way is to modify the
MR sequence by adding

#ifndef VXWORKS
	pMrProt->fwrite("E:\\image\\mrprot.txt")
#endif

to the function fSeqCheck. This will dump the complete protocol
information to the specified location just before the first scan.
With the PixelDataGrabber listening for files in E:\image, it will
note this and immediately parse the new protocol. The filename
"mrprot.txt" is currently hard-coded in the PixelDataGrabber.

If you run an unmodified MR sequence that does not dump the information,
you can try to create your own mrprot.txt and place it in e:\image before
running the scans. In this case, the PixelDataGrabber will read that file
when the first PixelData arrives. See section 3 for more info on the
protocol file.

In case the PixelDataGrabber encounters a mismatch between protocol
specifications and the size of the .PixelData files, it will report
an error and not write the sample.


3) Parsing protocol information using "siemensap"
-------------------------------------------------
"Siemensap" provides some plain C functions and datatypes to parse
the ASCII format Siemens protocol data into a list or tree of
key/value items. Currently supported value types are strings, long
integers, and double precision numbers. Field types are determined
automatically to a large extend (e.g. a dot implies a double), but
some special rules are added. For example, a field name that starts
with "d" will always be parsed as a double precision value.

If you need to write your own protocol file (E:\image\mrprot.txt),
the most important ingredients are the following lines:
alTR = 2900000
lContrasts = 5
sKSpace.lBaseResolution = 64
sSliceArray.lSize = 32
sSliceArray.asSlice[0].dPhaseFOV = 224.0
sSliceArray.asSlice[0].dReadoutFOV = 168.0
sSliceArray.asSlice[0].dThickness = 3.0
which are accessed in sap_get_essentials(...).

Please see siemensap.h for Doxygen-style documentation of the API.


4) fMRI visualisation clients
-----------------------------
See gui_buffer_client and opengl_client, both based on FLTK, 
platform-independent, work in progress. Might be handy for demos
and online inspection.


5) Integration in Matlab
------------------------
A demo is available: realtime_fmriviewer.m
Otherwise see usual FieldTrip documentation, e.g., ft_read_data and ft_read_header.
