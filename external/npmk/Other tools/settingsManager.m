function varargout = settingsManager(varargin)

% settingsManager
%
% Loads, maintains and saves the NPMK Settings File. For internal use only.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Use settingsManager
%
% Use settingsFile = settingsManager(settingsFile)
%
%   settingsFile:   The settings variable passed to the function for
%                   saving. If no input has been specified, the function
%                   will load an existing file and pass it out to the
%                   calling function. If a setting file does not exist, the
%                   function will generate a file with defalt settings and
%                   then will pass the variables to the calling function. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   support@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0: August 13, 2020
%   - Initial release.
%
% 1.1.0.0: September 16, 2020
%   - Fixed a bug where the response to 1/4 uV warning was not being saved.
%

% Finding the correct settings file for this version of MATLAB
fullPath = which('settingsManager.m');
fullPath = [fullPath(1:end-17) 'settings.npmk'];

% Does the file already exist?
if exist(fullPath, 'file') == 2
    if length(varargin) == 0 % Load the settings file and send as output
        varargout{1} = load(fullPath, '-mat');
        try varargout{1}.checkeddate; catch varargout{1}.checkeddate = datetime; end
        try varargout{1}.ShowZeroPadWarning; catch varargout{1}.ShowZeroPadWarning = 1; end
        try varargout{1}.ShowuVWarning; catch varargout{1}.ShowuVWarning = 1; end
    elseif length(varargin) == 1
        checkeddate = varargin{1}.checkeddate;
        ShowZeroPadWarning = varargin{1}.ShowZeroPadWarning;
        ShowuVWarning = varargin{1}.ShowuVWarning;
        save(fullPath, 'checkeddate', 'ShowZeroPadWarning', 'ShowuVWarning'); % Save the settings file
    end
else % Since it doesn't exist, create and save it.
    checkeddate = datetime;
    ShowZeroPadWarning = 1;
    save(fullPath, 'checkeddate', 'ShowZeroPadWarning');
    varargout{1} = load(fullPath, '-mat');
end


end % END of Function
