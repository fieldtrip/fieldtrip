function [data] = mne_ex_average_epochs(dataname,origname,outname)
%
%   An example of averaging over epochs
%
%   function mne_ex_average_epochs(dataname,origname,outname)
%
%   dataname  - Name of a epoch data. The description file will
%               be <dataname>_desc.mat and the epoch file <dataname>.epochs
%   origname  - Name of the file from which the epochs were extracted.
%   outname   - Name of the output file (optional)
%
%   Returns an evoked data structure identical to the ones 
%   returned from fiff_read_evoked
%   
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.2  2007/11/13 10:55:32  msh
%   Specify the second argument to all calls to the exist function
%
%   Revision 1.1  2006/09/27 20:15:11  msh
%   Added a new example to demonstrate averaging.
%
%

me='MNE:mne_ex_average_epochs';

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end

if nargin ~= 2 && nargin ~= 3
    error(me,'Incorrect number of arguments');
end
%
%   Set up file names
%
descname  = sprintf('%s_desc.mat',dataname);
epochname = sprintf('%s.epochs',dataname);
%
%   Check that all the necessary files exist
%
if exist(descname,'file') ~= 2
    error(me,'The description file %s does not exist',descname);
end
if exist(epochname,'file') ~= 2
    error(me,'The epoch data file %s does not exist',epochname);
end
if exist(origname,'file') ~= 2
    error(me,'The original data file %s does not exist',origname);
end
%
%   Load the epoch info
%
load(descname)
fprintf(1,'Epoch info loaded\n');
%
%   Read the measurement info from the original data file
%
try
    info = fiff_read_meas_info(origname);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end
fprintf(1,'Measurement info loaded\n');
%
%   Adjust the measurement info with the information from epoch description
%
for k = 1:MNE_epoch_info.nchan
    sel{k} = deblank(MNE_epoch_info.ch_names(k,:));
end
info = fiff_pick_info(info,fiff_pick_channels(info.ch_names,sel));

info.sfreq = MNE_epoch_info.sfreq;
info.lowpass = MNE_epoch_info.lowpass;
info.highpass = MNE_epoch_info.highpass;

fprintf(1,'Averaging');
%
%   Load the first epoch
%
try
    [one,fid] = mne_read_epoch(MNE_epoch_info,1,-1);
catch
    error(me,'Failed to read the first epoch');
end
ave = one;
nave = 1;
fprintf(1,'.');
%
%   Load all subsequent ones and average
%
for k = 2:MNE_epoch_info.nepoch
    try
        [one,fid] = mne_read_epoch(MNE_epoch_info,k,fid);
    catch
        error(me,'Failed to read the %dth epoch',k);
    end
    %
    %   You should design an artefact rejection routine yourself   
    %
    ave = ave + one;
    nave = nave + 1;
    fprintf(1,'.');
end
ave = ave/nave;
fprintf(1,'[done]\n');
%
%   Just put together...
%
evoked.aspect_kind = int32(FIFF.FIFFV_ASPECT_AVERAGE);
evoked.nave        = int32(nave);
%
%   This makes the reasonable assumption that all epochs are 
%   not only of the same length but also have the same time scale
%
evoked.first   = int32(MNE_epoch_info.epochs(1,4));
evoked.last    = int32(evoked.first + MNE_epoch_info.epochs(1,5) - 1);
evoked.comment = 'Put your own comment here';
evoked.times   = double(evoked.first:1:evoked.last)/info.sfreq;
evoked.epochs  = ave;
%
data.info   = info;
data.evoked = evoked;
%
%   Do we want the output?
%
if nargin == 3
    try 
        fiff_write_evoked(outname,data);
        fprintf(1,'Wrote %s\n',outname);
    catch
        error(me,'Failed to write the data');
    end
end


    

