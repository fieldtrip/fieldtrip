function NSxToHL(fname)

% NSxToHL
% 
% Opens and reads an NSx file without the header information and saves
% the binary data into a .dat file with the same name. This can be used for
% specific applications that require this type of data, e.g. Klusters. 
% Works with File Spec 2.1, 2.2 and 2.3.
%
% It does not support pauses at this time.
%
% Use OUTPUT = NSxToHL(fname)
% 
% All input arguments are optional.
%
%   fname:        Name of the file to be opened. If the fname is omitted
%                 the user will be prompted to select a file. 
%                 DEFAULT: Will open Open File UI.
%
%   OUTPUT:       Contains the binary data.
%
%   Example 1: 
%   NSxToHL('c:\data\sample.ns5');
%
%   In the example above, the file c:\data\sample.ns5 will be opened and
%   the data will be read and saved without the header information in file
%   sample.ns5.dat.
%
%   Example 2:
%   NSxToHL;
%
%   In the example above, the file user will be prompted for the file. The
%   selected file (e.g. FILENAME.NSx) will be opened and the data will be
%   read and saved in the file FILENAME.NSx.dat.
%
%   Kian Torab
%   kian@blackrockmicro.com
%   Blackrock Microsystems
%   Version 1.0.0.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0:
%   - Initial release.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Opening the file
% Popup the Open File UI. Also, process the file name, path, and extension
% for later use, and validate the entry.
if ~exist('fname', 'var')
    if ~ismac
        [fname, path] = getFile('*.ns*', 'Choose an NSx file...');
    else
        [fname, path] = getFile('*.*', 'Choose an NSx file...');
    end
    if fname == 0
        disp('No file was selected.');
        if nargout
            clear variables;
        end
        return;
    end
else
    if isempty(fileparts(fname))
        fname = which(fname);
    end
    [path,fname, fext] = fileparts(fname);
    fname = [fname fext];
    path  = [path '/'];
end
if fname==0
    return; 
end

%% Reading the headerless file
data = openNSxHL([path fname]);

%% Writing to file
% Determining the filename for the converted file
newFilename = [path fname '.dat'];

% Opening the output file for saving
FIDw = fopen(newFilename, 'w+', 'ieee-le');

% Writing data into file
disp('Writing the converted data into the new .dat file...');
fwrite(FIDw, data, 'int16');