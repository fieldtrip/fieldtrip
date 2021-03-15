# Known limitations to the EGI/Philips MFF import/export plugin
---------------------------------------------------------
- Timing accuracy limited to 10-50 microseconds precision. Matlab 2014b time 
  conversion is not as precise as Matlab 2017a conversion (we can loose 1ms 
  accuracy on some events in some specific rare files). With Matlab 2018a, the
  difference seems to be at most 10 microseconds. In most applications,
  this limitation has no consequences as event latencies are multiple of 1000
  microseconds (1 ms) and the events are imported and exported perfectly. 
 
- The plugin will not work with version of Matlab older than 2014a as the Java 
  JAR file cannot be properly interfaced.

- When the Java Heap Memory is at its default level, few files can be imported. 
  There is a special message guiding users how to increase the default settings 
  in Matlab in case they encounter an error. Importantly, some large files
  require a computer with 16Gb of RAM or more (or Java Heap of 4Gb).

- Calibration data when present is applied to the data. However, if the file
  is exported this information is lost. Impedance data and channel status also 
  contained in the info1.xml file is ignored.

- Importing trials of different length is not supported. This is a rare occurrence.

- Channel status and channel keys in category files are not imported

- Importing filter information and calibration information is not supported
  (except gain calibration which is applied to the data when it is imported)

- The plugin will import multiple subject trial averages in EEGLAB but the 
   export is not correct or identical. The reason is that EEGLAB consider 
   these averages as single trials from the same subject and cannot handle 
   participant names.

- When importing in standalone mode (not as an EEGLAB plugin, if the file 
  is exported, it will not contain event keys

- Video files are not imported and expoted

- This plugin was tested on platforms using little-endian byte ordering. 
  Although we do not expect big-endian to be a problem, there is a small
  chance there could be problem

# Revision history
----------------
Version 3.4
- Fixes issues associated with reading data under Fieldtrip

Version 3.3
- Fix reading on Windows

Version 3.2
- Fix error with boundary events

Version 3.1
- Implement auto-save for BIDS
- Fix event track type for segmenting in Netstation

Version 3.0
- Fix issue with PNS channels
- Fix imporing rare files with custom fields
- Fixing a variety of minor Git issues

Version 2.02
- Adding eeg_checkset to private folder to run in standalone mode (F. Tadel)
- Added fix for files with no events (R. Martinez fix for https://github.com/sccn/eeglab/issues/17)

Version 2.01
- Allow eegplugin_mffmatlabio to return version number
- Remove call in mff_import that was assuming EEGLAB was present

Version 2.00
- Octave compatibility
- Fix issue with boundary latency when importing file mff version 0
- Allowing to export random EEG files
- Rescale coordinates for non-MFF channel coordinates
- Allow exporting datasets which do not have a code field
- Allow exporting datasets with no event duration
- Better support for PNS channels for File-IO
- Fix command line call not rotating channels

Version 1.00
- Add file separator to EEGLAB export menu

Version 0.96
- Add license for each file
- Clean up documentation

Version 0.95
- Fix issues when running File-io import and now importing using File-io functions direclty
- Adding licence file

Version 0.94
- Fix EEGLAB history for pop_mffimport
- Fix boundaries when encoding types
- Minor documentation changes

Version 0.93
- Now import/export data files with PNS data only
- Through the graphic interface, now allow to specify the MFF events field to
  use for the EEGLAB event types.
- Now allow the plugin to function in standalone mode.
- Now check for Matlab version and issue an error for unsuported Matlab versions

Version 0.92
- Renamed all the functions
- Fix minor issue with multiple references

Version 0.91 - Difference with previous revision
- Fixed issue with info1.xml file which was missing some information

Version 0.9
- Adding support for layout and subject information

