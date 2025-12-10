function timeDelay = calcTimeDelay(NSx1, NSx2)

%
% timeDelay = calcTimeDelay(NSx1, NSx2)
%
%   This function calculates the time difference between the end of the
%   first NSx file and the beginning of the second NSx file.
%
%   NSx1:       The first NSx file passed to the function. This input is 
%               optional. In its absense, a dialog will open and will 
%               prompt the user to select an NSx file.
%               (OPTIONAL)
%
%   NSx2:       The second NSx file passed to the function. This input is 
%               also optional. In its absense, a dialog will open and will
%               prompt the user to select an NSx file.
%               (OPTIONAL)
%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%
%   Version 1.0.0.0

if ~exist('NSx1', 'var') || ~exist('NSx2', 'var')
    NSx1 = openNSx;
    NSx2 = openNSx;
end

timeDifference = double(NSx2.MetaTags.DateTimeRaw - NSx1.MetaTags.DateTimeRaw);
timeDifference = timeDifference(end-1)+timeDifference(end)/1000;

timeDelay = (timeDifference*NSx1.MetaTags.SamplingFreq-NSx1.MetaTags.DataPoints)/NSx1.MetaTags.SamplingFreq;

if ~nargout
    disp(['The time delay between the end of the first file and the beginning of the second file is ' num2str(timeDelay) ' seconds.']);
    return;
end