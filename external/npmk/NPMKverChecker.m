function NPMKverChecker()

% NPMKverChecker
%
% Checks to see if there is a newer version of NPMK is available for
% download.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Use NPMKverChecker
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   support@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0: September 13, 2017
%   - Initial release.
%
% 1.0.1.0: September 13, 2017
%   - Fixed a crash in case there is no Internet connection.
%
% 1.0.2.0: January 10, 2018
%   - Added a clickable URL to the prompt.
%
% 1.1.0.0: January 27, 2020
%   - Only checks for a new version once a week instead of every time.
%

%% Variables and constants
gitHubURL = 'https://github.com/BlackrockMicrosystems/NPMK/releases/latest';

%% Find full path of NPMKverChecker.m
fileFullPath = which('NPMKverChecker.m');
fileFullPath = [fileFullPath(1:end-1) 'dat'];

%% Check for the latest version fo NPMK
try
    if exist(fileFullPath, 'file') == 2
        load(fileFullPath, '-mat');
        if floor(abs(now - datenum(checkeddate - days(1)))) > 8 %#ok<NODEF>
            disp('Checking for a new version of NPMK...');
            checkver = 1;
        else
            checkver = 0;
        end
    else
        checkver = 1;
    end
    if checkver
        FIDv = fopen('Versions.txt');
        verFile = fscanf(FIDv, '%s'); 
        fclose(FIDv);
        latestVersion = verFile(findstr('LATEST', verFile)+7:findstr('LATEST', verFile)+13);
        gitHubPage = urlread(gitHubURL);
        newVersionAvailable = findstr(latestVersion, gitHubPage);
        if isempty(newVersionAvailable)
            disp('A new version of NPMK may be available.');
            fprintf('Please visit <a href="%s">GitHub NPMK Page</a> to get the latest version.\n', gitHubURL)
        end
        checkeddate = datetime;
        save(fileFullPath, 'checkeddate');
    end
catch
end