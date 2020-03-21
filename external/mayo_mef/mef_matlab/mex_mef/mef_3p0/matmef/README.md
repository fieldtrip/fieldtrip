# matmef
Matlab wrapper for MEF 3.0 library

## Introduction
Several Matlab mex functions that wrap around the MEF 3.0 library to read MEF 3.0 data

- Written by Max van den Boom (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)
- Adapted from PyMef (by Jan Cimbalnik, Matt Stead, Ben Brinkmann, and Dan Crepeau)

## Using pre-compiled mex files
1. Clone or download (and extract) the matmef repository
2. Use the functions

## Building from source
1. Clone the matmef repository using: git clone --recurse-submodules https://github.com/MaxvandenBoom/matmef.git
   - Note: unfortunately git clone does not clone submodules by default, so make sure to add '--recurse-submodules' option
2. Start matlab and set the matmef folder as your working directory
3. To compile the .mex files, run the following lines in matlab:
   - mex read_mef_session_metadata.c meflib/meflib/meflib.c meflib/meflib/mefrec.c matmef_mapping.c mex_datahelper.c
   - mex read_mef_ts_data.c matmef_data.c meflib/meflib/meflib.c meflib/meflib/mefrec.c

## Matlab usage examples
%  
% High-level  
%  

% metadata only  
[metadata] = readMef3('./mefSessionData/');  

% two channels  
[metadata, signaldata] = readMef3('./mefSessionData/', [], {'Ch02', 'Ch07'});  

% all channels, samples 0-1000  
[metadata, signaldata] = readMef3('./mefSessionData/', [], [], 'samples', 0, 1000);  

% two channels, samples 0-1000  
[metadata, signaldata] = readMef3('./mefSessionData/', [], {'Ch02', 'Ch07'}, 'samples', 0, 1000);  
  
  
%  
% Low-level  
%  
  
session = read_mef_session_metadata('./mefSessionData/', [], 1);  
data = read_mef_ts_data('./mefSessionData/channelPath/');  
data = read_mef_ts_data('./mefSessionData/channelPath/', [], 'samples', 0, 1000);  
data = read_mef_ts_data('./mefSessionData/channelPath/', [], 'time', 1578715810000000, 1578715832000000);  
