function NSx = openNSxSync(varargin)

% openNSxSync
% 
% Opens a synced NSx file and removed the extra bit of data from the file.
%
% INPUT
%
%   filename:   The name of the file to be opened.
%               DEFAULT: The program will prompt the user to select a file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   kianblackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0:
%   - Initial release.
%
% 1.1.0.0: June 13, 2014
%   - Added the ability to open a file by passing on the file name.
%
% 1.1.1.0: December 3, 2014
%   - Fixed a minor bug.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializing varilables.
for inputIDX = 1:nargin
    inputVar = varargin{inputIDX};
    if isnumeric(inputVar)
        if inputVar >= 1 && inputVar <= 2
            method = inputVar;
        else
            disp('The variable method can only be 1 (defauly) or 2. See help for more information.');
            NSx = -1;
            return;
        end
    else
        filename = inputVar;
    end
end

%% Openning Synced files and removing the extra piece of data
if exist('filename', 'var')
    if exist(filename, 'file') == 2
        NSx = openNSx(filename);
    else
        disp('File was not found.');
        NSx = -1;
        return;
    end
else
    NSx = openNSx;
end

%% Getting rid of the data extra bits and pieces created by the sync algorithm
if iscell(NSx.Data)
    % Removing the extra bit of empty data
    NSx.Data = NSx.Data{2};
    NSx.MetaTags.Timestamp(1) = [];
    NSx.MetaTags.DataPoints(1) = [];
    NSx.MetaTags.DataDurationSec(1) = [];
    % Re-aligning what's left
    NSx.Data = [zeros(NSx.MetaTags.ChannelCount, floor(NSx.MetaTags.Timestamp)/NSx.MetaTags.SamplingFreq) NSx.Data];
    NSx.MetaTags.DataPoints = NSx.MetaTags.DataPoints - NSx.MetaTags.Timestamp;
    NSx.MetaTags.DataDurationSec = NSx.MetaTags.DataPoints / NSx.MetaTags.SamplingFreq;
    NSx.MetaTags.Timestamp = 0;
end

%% If user does not specify an output argument it will automatically create a structure.
outputName = ['NS' NSx.MetaTags.FileExt(4)];
if (nargout == 0),
    assignin('caller', outputName, NSx);
    clear all;
else
    varargout{1} = NSx;
end
