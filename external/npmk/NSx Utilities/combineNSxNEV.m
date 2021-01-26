%
%   [NSx, NEV] = combineNSxNEV(filename1, filename2)
%
%   This function loads two NSx and NEV files and it will combine them
%   together into a sinlge NSx and NEV structure in MATLAB. To merge two
%   NSx and NEV files into individual files see mergeNSxNEV. The time
%   difference between the two sets of recordings is removed. To determine
%   the time differnce between the two data files, use
%   NSx.MetaTags.DateTimeRaw or NEV.MetaTags.DateTimeRaw variables.
%
%
%   filename1:  The name of the first NSx file. This input is optional. In
%               its absense, a dialog will open and will prompt the user to
%               select an NSx file.
%               (OPTIONAL)
%
%   filename2:  The name of the second NSx file. This input is also
%               optional. In its absense, a dialog will open and will
%               prompt the user to select an NSx file.
%               (OPTIONAL)
%
%
%   Example: 
%   
%   [NSx, NEV] = combineNSxNEV('c:\data\saveddata1.ns5', 'c:\data\saveddata2.ns5');
%
%   The above example reads the two files (full path needed)
%   c:\data\saveddata1.ns5 and c:\data\saveddata2.ns5 and their corresponding
%   NEV files (saveddata1.nev and saveddata2.nev) in the same folder and
%   combines them into single variables NSx and NEV into MATLAB workspace.
%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%
%   Version 1.1.2.0


function [NSx1, NEV1] = combineNSxNEV(filename1, filename2)

% Openning NSx files
if exist('filename1', 'var') && exist('filename2', 'var')
    disp('Load the first NSx file.');
    NSx1 = openNSx('read', filename1);
    disp('Load the second NSx file.');
    NSx2 = openNSx('read', filename2);
    if NSx1.MetaTags.SamplingFreq ~= NSx2.MetaTags.SamplingFreq;
        disp('The sampling frequencies are not the same.');
        return;
    end
else
    NSx1 = openNSx('read');
    NSx2 = openNSx('read');
end

% Determining length of the first NSx file
conversionFactor = 30000/NSx1.MetaTags.SamplingFreq;
NSx1DataLength = NSx1.MetaTags.DataPoints * conversionFactor;

% Combining NSx files
NSx1.Data = [NSx1.Data, NSx2.Data];

% Opening NEV files
fileNameNEV1 = [NSx1.MetaTags.FilePath '/' NSx1.MetaTags.Filename(1:end-3), 'nev'];
fileNameNEV2 = [NSx1.MetaTags.FilePath '/' NSx2.MetaTags.Filename(1:end-3), 'nev'];
clear NSx2;

if (exist(fileNameNEV1, 'file') == 2) && (exist(fileNameNEV2, 'file') ==2)
    disp('Openning corresponding NEV files...');
    NEV1 = openNEV('read', fileNameNEV1);
    NEV2 = openNEV('read', fileNameNEV2);
else
    disp('Load the first NEV file.');
    NEV1 = openNEV('read');
    disp('Load the second NEV file.');    
    NEV2 = openNEV('read');
end

% Adjusting the timestamp on the second NEV file
NEV2.Data.Comments.TimeStamp = NEV2.Data.Comments.TimeStamp + NSx1DataLength;
NEV2.Data.Comments.TimeStampSec = NEV2.Data.Comments.TimeStampSec + double(NSx1DataLength)/30;
NEV2.Data.SerialDigitalIO.TimeStamp = NEV2.Data.SerialDigitalIO.TimeStamp + NSx1DataLength;
NEV2.Data.SerialDigitalIO.TimeStampSec = NEV2.Data.SerialDigitalIO.TimeStampSec + double(NSx1DataLength)/30;
NEV2.Data.Spikes.TimeStamp = NEV2.Data.Spikes.TimeStamp + NSx1DataLength;
NEV2.Data.VideoSync.TimeStamp = NEV2.Data.VideoSync.TimeStamp + NSx1DataLength;
if ~isempty(NEV2.Data.Tracking)
    trackingFieldNames = fieldnames(NEV2.Data.Tracking);
    for idx = 1:size(trackingFieldNames, 1)
        NEV2.Data.Tracking.(trackingFieldNames{idx}).TimeStamp = NEV2.Data.Tracking.(trackingFieldNames{idx}).TimeStamp + NSx1DataLength;
    end
end
NEV2.Data.PatientTrigger.TimeStamp = NEV2.Data.PatientTrigger.TimeStamp + NSx1DataLength;
NEV2.Data.Reconfig.TimeStamp = NEV2.Data.Reconfig.TimeStamp + NSx1DataLength;

% Combining the two NEV files
NEV1.Data.Spikes.Electrode      = [NEV1.Data.Spikes.Electrode, NEV2.Data.Spikes.Electrode];
NEV1.Data.Spikes.TimeStamp      = [NEV1.Data.Spikes.TimeStamp, NEV2.Data.Spikes.TimeStamp];
NEV1.Data.Spikes.Unit           = [NEV1.Data.Spikes.Unit, NEV2.Data.Spikes.Unit];
NEV1.Data.Spikes.Waveform       = [NEV1.Data.Spikes.Waveform, NEV2.Data.Spikes.Waveform];
NEV1.Data.Comments.TimeStamp    = [NEV1.Data.Comments.TimeStamp, NEV2.Data.Comments.TimeStamp];
NEV1.Data.Comments.TimeStampSec = [NEV1.Data.Comments.TimeStampSec, NEV2.Data.Comments.TimeStampSec];
NEV1.Data.Comments.CharSet      = [NEV1.Data.Comments.CharSet, NEV2.Data.Comments.CharSet];
NEV1.Data.Comments.Color        = [NEV1.Data.Comments.Color, NEV2.Data.Comments.Color];
NEV1.Data.Comments.Text         = [NEV1.Data.Comments.Text; NEV2.Data.Comments.Text];
if ~isempty(NEV2.Data.Tracking)
    for idx = 1:size(trackingFieldNames, 1)
        NEV1.Data.Tracking.(trackingFieldNames{idx}).TimeStamp = [NEV1.Data.Tracking.(trackingFieldNames{idx}).TimeStamp, NEV2.Data.Tracking.(trackingFieldNames{idx}).TimeStamp];
        NEV1.Data.Tracking.(trackingFieldNames{idx}).TimeStampSec = [NEV1.Data.Tracking.(trackingFieldNames{idx}).TimeStampSec, NEV2.Data.Tracking.(trackingFieldNames{idx}).TimeStampSec];
        NEV1.Data.Tracking.(trackingFieldNames{idx}).ParentID = [NEV1.Data.Tracking.(trackingFieldNames{idx}).ParentID, NEV2.Data.Tracking.(trackingFieldNames{idx}).ParentID];
        NEV1.Data.Tracking.(trackingFieldNames{idx}).NodeCount = [NEV1.Data.Tracking.(trackingFieldNames{idx}).NodeCount, NEV2.Data.Tracking.(trackingFieldNames{idx}).NodeCount];
        NEV1.Data.Tracking.(trackingFieldNames{idx}).MarkerCount = [NEV1.Data.Tracking.(trackingFieldNames{idx}).MarkerCount, NEV2.Data.Tracking.(trackingFieldNames{idx}).MarkerCount];
        NEV1.Data.Tracking.(trackingFieldNames{idx}).MarkerCoordinates = [NEV1.Data.Tracking.(trackingFieldNames{idx}).MarkerCoordinates NEV2.Data.Tracking.(trackingFieldNames{idx}).MarkerCoordinates];
    end
end