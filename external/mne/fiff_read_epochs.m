function [epochs] = fiff_read_epochs(fname)
%
% [epochs] = fiff_read_epochs(fname,setno)
%
% Read eochs from file
%
%
%   Author : Martin Luessi, MGH Martinos Center
%   License : BSD 3-clause

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me='MNE:fiff_read_epochs';

%
%   Open the file
%
fprintf(1,'Reading %s ...\n',fname);
[ fid, tree ] = fiff_open(fname);
%
%   Read the measurement info
%
[ info, meas ] = fiff_read_meas_info(fid,tree);
info.filename = fname;

%
%   Read events
%
[ events, mappings ] = fiff_read_events(fid,tree);

%
% Locate the data of interest
%
processed = fiff_dir_tree_find(meas, FIFF.FIFFB_PROCESSED_DATA);
if length(processed) == 0
        fclose(fid);
    error(me,'Could not find epochs data');
end

ep = fiff_dir_tree_find(meas, FIFF.FIFFB_MNE_EPOCHS);
if length(ep) == 0
        fclose(fid);
    error(me,'Could not find epochs data');
end

comment = '';
selection = '';
drop_log = '';
for k = 1:ep.nent
    kind = ep.dir(k).kind;
    pos  = ep.dir(k).pos;
    switch kind
        case FIFF.FIFF_FIRST_SAMPLE
            tag = fiff_read_tag(fid,pos);
            first = tag.data;
        case FIFF.FIFF_LAST_SAMPLE
            tag = fiff_read_tag(fid,pos);
            last = tag.data;
        case FIFF.FIFF_COMMENT
            tag = fiff_read_tag(fid,pos);
            comment = tag.data;
        case FIFF.FIFF_EPOCH
            tag = fiff_read_tag(fid,pos);
            epoch = tag.data;
        case FIFF.FIFF_MNE_BASELINE_MIN
            tag = fiff_read_tag(fid,pos);
            bmin = tag.data;
        case FIFF.FIFF_MNE_BASELINE_MAX
            tag = fiff_read_tag(fid,pos);
            bmax = tag.data;
        case FIFF.FIFFB_MNE_EPOCHS_SELECTION
            tag = fiff_read_tag(fid,pos);
            selection = tag.data;
        case FIFF.FIFFB_MNE_EPOCHS_DROP_LOG
            tag = fiff_read_tag(fid,pos);
            drop_log = tag.data;
    end
end

if ~exist('epoch','var')
    fclose(fid);
    error(me,'Epochs data not found');
end

if ~exist('bmin','var')
    bmin = double(first)/info.sfreq;
end

if ~exist('bmax','var')
    bmax = double(last)/info.sfreq;
end

baseline = [bmin, bmax];

nsamp = last-first+1;
fprintf(1,'\tFound the data of interest:\n');
fprintf(1,'\t\tt = %10.2f ... %10.2f ms (%s)\n',...
    1000*double(first)/info.sfreq,1000*double(last)/info.sfreq,comment);
if ~isempty(info.comps)
    fprintf(1,'\t\t%d CTF compensation matrices available\n',length(info.comps));
end

if size(epoch,1) ~= size(events,1)
    fclose(fid);
    error(me,'Incorrect number of trials (%d instead of %d)',...
        size(epoch,1),size(events,1));
end

if size(epoch,2) ~= info.nchan
    fclose(fid);
    error(me,'Incorrect number of channels (%d instead of %d)',...
        size(epoch,2),info.nchan);
end

if size(epoch,3) ~= nsamp
    fclose(fid);
    error(me,'Incorrect number of samples (%d instead of %d)',...
        size(epoch,3),nsamp);
end

%
%   Calibrate
%
for k = 1:info.nchan
    cals(k) = info.chs(k).cal;
end

nepochs = size(epoch, 1);
epoch = repmat(cals, [nepochs, 1, nsamp]) .* epoch;

times = double(first:last) / info.sfreq;
tmin = times(1);
tmax = times(end);

%
% Put it all together
%
epochs.info = info;
epochs.events = events;
epochs.name = comment;
epochs.times = times;
epochs.tmin = tmin;
epochs.tmax = tmax;
epochs.data = epoch;
epochs.baseline = baseline;
epochs.event_id = mappings;
epochs.selection = selection;
epochs.drop_log = drop_log;

fclose(fid);

return;
