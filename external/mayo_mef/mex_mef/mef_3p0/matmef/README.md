# matmef
Matlab wrapper for MEF 3.0 library

This project provides several Matlab mex functions that wrap around the MEF 3.0 library to read MEF 3.0 data

## Using pre-compiled mex files
1. Clone or download (and extract) the matmef repository
2. Use the functions

## Building from source
1. Clone the matmef repository using: `git clone https://github.com/MaxvandenBoom/matmef.git`
2. Start matlab and set the matmef folder as your working directory
3. To compile the .mex files, run the following lines in matlab:

   - `mex read_mef_session_metadata.c matmef_mapping.c mex_datahelper.c`
   - `mex read_mef_ts_data.c matmef_data.c`

## Matlab usage examples
```
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

% all channels, multiple ranges/epochs
ranges = [[0,    1000]; ...
          [1000, 2000]; ...
          [5000, 6000]];
[metadata, signaldata] = readMef3('./mefSessionData/', [], [], 'samples', ranges);
```


```
%  
% Low-level  
%  
  
session = read_mef_session_metadata('./mefSessionData/', [], 1);  
data = read_mef_ts_data('./mefSessionData/channelPath/');  
data = read_mef_ts_data('./mefSessionData/channelPath/', [], 'samples', 0, 1000);  
data = read_mef_ts_data('./mefSessionData/channelPath/', [], 'time', 1578715810000000, 1578715832000000);  
```

## Acknowledgements

- Written by Max van den Boom (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)
- Adapted from PyMef (by Jan Cimbalnik, Matt Stead, Ben Brinkmann, and Dan Crepeau)

- This project was funded by the National Institute Of Mental Health of the National Institutes of Health Award Number R01MH122258 to Dora Hermes
