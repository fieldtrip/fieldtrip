function codeKeeper = syncPatternFinderNSx(filenameNSx)

% syncPatternFinderNSx
% 
% Canculated the unique 8-bit code and the associated timestamp of the SYNC
% pulse of the NSP. 
%
% Use codeKeeper = syncPatternFinderNSx(filenameNSx)
% 
% All input arguments are optional. Input arguments can be in any order.
%
%   filenameNSx:  Name of the file to be opened. If the fname is omitted
%                 the user will be prompted to select a file. 
%                 DEFAULT: Will open Open File UI.
%
%
%   OUTPUT:       The structure that contains all the unique 8-bit SYNC
%                 pulse codes and their corresponding timestamps.
%
%   Example 1: 
%   codeKeeper = syncPatternFinderNSx('myTestNS5.ns5');
%
%   In the example above, the file myTestNS5.ns5 will be opened and all the
%   SYNC pulses and their corresponding timestamps will be output in the
%   codeKeeper structure. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   kian@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version History
%
% 1.0.0.0:  November 1, 2014
%   - Initial release.
%
% 1.0.1.0:  March 28, 2015
%   - Fixed a bug where the extension wsan't specified for the file when
%     the input file was not specified.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Opening the file to get the header information
if ~exist('filenameNSx', 'var')
    NSx = openNSx('noread');
    % Determining the file name read
    filenameNSx = fullfile(NSx.MetaTags.FilePath, NSx.MetaTags.Filename, NSx.MetaTags.FileExt);
else
    NSx = openNSx(filenameNSx, 'noread');
end

if ~isstruct(NSx)
    disp('Error reading the NSx file. Most likely the file name did not exist.');
    return;
end

%% Calculating the variables for the sync signal
% Figuring out the channel that recorded the sync pulse. This assums that
% the sync pulse was recorded on analog input 16 on the Cerebus or analog
% input 3 on the Direct system.
numberOfChannels = NSx.MetaTags.ChannelCount;

% Calculating the maximum number of points to read
maxPacketsToRead = 30 * NSx.MetaTags.SamplingFreq;
if maxPacketsToRead > NSx.MetaTags.DataPoints
    maxPacketsToRead = NSx.MetaTags.DataPoints;
end

%% Reading the sync signal
NSxSync = openNSx(filenameNSx, ['c:' num2str(numberOfChannels)], ['t:1:' num2str(maxPacketsToRead)], 'read');

codeKeeper = syncPatternDetectNSx(NSxSync.Data(1,:));