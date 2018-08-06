function datastruct = ghdf5read(filename)
%
%GHDF5READ reads data recorded with the g.Recorder data acquisition software.
%    DATASTRUCT = GHDF5READ(FILENAME) reads data in HDF5 meta file format.
%    For more information see http://www.hdfgroup.org. 
%    ghdf5read.m is based on ghdf5fileimport.p.
%
%    FILENAME: string containing the name of file to import.
%
%    DATASTRUCT: MATLAB structure containing all information stored in
%    the hdf5 file. This structure may contain various sub-structures, depending on
%    the specific recording settings.
%
%    Example:
%    Read in the raw 'data' and the used sampling frequency 'fs' from:
%    filename = 'RecordSession_2009.10.29_09.43.20.hdf5';
%    datastruct = ghdf5read(filename);
%    data = datastruct.RawData.Samples;
%    fs = datastruct.RawData.AcquisitionTaskDescription.SamplingFrequency;
%
%   (c) g.tec medical engineering GmbH
% 

% To access the raw data and device dependent information the following
% sub-structures are important:
%
% datastruct.AsynchronData
% datastruct.RawData
% datastruct.Version
%
% AsynchronData:
% --------------
% Trigger and Marker information are found in 
% AsynchronSignalTypes: Containing one sub-structure AsynchronSignalDescription
%                       for each marker or trigger with name, ID and other information of
%                       selected triggers and markers.
% Time:                 sample number of marker or trigger.
% TypeID:               ID corresponding to ID in AsynchronSignalDescription.
% Value:                high (1) or low (0).
%
% RawData:
% --------
% Raw data information and samples:
% AcquisitionTaskDescription: Contains data about the acquired channels
%                        such as sensitivity, sample rate, etc.This sub group
%                        describes the data stored in the Samples array.
% DAQDeviceCapabilities: Contains hardware specific data about all channels of the device.
% DAQDeviceDescription:  Some information about the device such as device name,
%                        vendor name, …
% Samples:               Contains the data of all aquired channels.
% SessionDescription:    Describes the recording session.
% SubjectDescription:    Describes the subject.
%
% Version:
% --------
% file version of hdf5 file.

datastruct = ghdf5fileimport(filename);