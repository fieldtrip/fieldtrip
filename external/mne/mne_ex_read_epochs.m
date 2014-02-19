function [data,times,ch_names] = mne_ex_read_epochs(fname,event,eventname,tmin,tmax)
%
%   Example of reading raw data
%
%   [ data, times, ch_names ] = mne_ex_read_epochs(fname,event,eventname,tmin,tmax)
%
%   Input :
%
%   fname       - The name of the input file
%   event       - The event
%   eventname   - Name of the event file
%   tmin        - Starting time in seconds
%   tmax        - Ending time in seconds
%
%   Output :
%
%   data        - Array of structures corresponding to the epochs with fields:
%
%                 epoch    the epoch, channel by channel
%                 event    event #
%                 tmin     starting time in the raw data file (initial skip omitted)
%                 tmax     ending stime in the raw data file (initial skip omitted)
%
%   times       - The time points of the samples, in seconds
%   ch_names    - Names of the channels included
%
%
%   NOTE 1: The purpose of this function is to demonstrate the raw data reading
%   routines. You may need to modify this for your purposes
%
%   NOTE 2: You need to run mne_process_raw once as
%
%   mne_process_raw --raw <fname> --projoff
%
%   to create the fif-format event file (or open the file in mne_browse_raw).
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.7  2009/03/21 00:17:35  msh
%   Added one missing semicolon.
%
%   Revision 1.6  2009/03/12 17:59:06  msh
%   Included channel names were not picked.
%
%   Revision 1.5  2009/01/27 22:02:17  msh
%   Added some more diagnostic output
%
%   Revision 1.4  2009/01/27 21:55:44  msh
%   Fixed error in reading a text event file with negative sample numbers.
%   The conversion from time to event numbers was not correctly done.
%
%   Revision 1.3  2009/01/12 20:30:16  msh
%   Improved the comments
%
%   Revision 1.2  2009/01/10 18:21:26  msh
%   Added an option to specify the event file
%   Changed the format of the output data
%
%   Revision 1.1  2008/10/10 16:13:57  msh
%   Added mne_ex_read_epochs. Fixed help text of mne_write_cov_file.m
%
%

%
%   Fiddle with the arguments
%
me='MNE:mne_ex_read_epochs';

keep_comp = false;
dest_comp = 0;
pick_all  = true;

if nargin ~= 5
    error(me,'Incorrect number of arguments');
end
%
%   Setup for reading the raw data
%
try
    raw = fiff_setup_read_raw(fname);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end

if pick_all
    %
    % Pick all
    %
    picks = 1:raw.info.nchan;
    %
else
    %
    %   Set up pick list: MEG + STI 014 - bad channels (modify to your needs)
    %
    include{1} = 'STI 014';
    want_meg   = true;
    want_eeg   = false;
    want_stim  = false;
    %
    %
    picks = fiff_pick_types(raw.info,want_meg,want_eeg,want_stim,include,raw.info.bads);
end
ch_names   = raw.info.ch_names(picks);
%
%   Set up projection
%
if isempty(raw.info.projs)
    fprintf(1,'No projector specified for these data\n');
    raw.proj = [];
else
    %
    %   Activate the projection items
    %
    for k = 1:length(raw.info.projs)
        raw.info.projs(k).active = true;
    end
    fprintf(1,'%d projection items activated\n',length(raw.info.projs));
    %
    %   Create the projector
    %
    [proj,nproj] = mne_make_projector_info(raw.info);
    if nproj == 0
        fprintf(1,'The projection vectors do not apply to these channels\n');
        raw.proj = [];
    else
        fprintf(1,'Created an SSP operator (subspace dimension = %d)\n',nproj);
        raw.proj = proj;
    end
end
%
%   Set up the CTF compensator
%
current_comp = mne_get_current_comp(raw.info);
if current_comp > 0
    fprintf(1,'Current compensation grade : %d\n',current_comp);
end
if keep_comp
    dest_comp = current_comp;
end
if current_comp ~= dest_comp
    try
        raw.comp = mne_make_compensator(raw.info,current_comp,dest_comp);
        fprintf(1,'Appropriate compensator added to change to grade %d.\n',dest_comp);
    catch
        error(me,'%s',mne_omit_first_line(lasterr));
    end
end
%
%  Read the events
%
if isempty(eventname)
    p = strfind(fname,'.fif');
    if p > 1
        eventname = sprintf('%s-eve.fif',fname(1:p-1));
    else
        error(me,'Raw file name does not end properly');
    end
    events = mne_read_events(eventname);
    fprintf(1,'Events read from %s\n',eventname);
else
    %
    %   Binary file
    %
    p = strfind(eventname,'-eve.fif');
    if p > 1
        try
            events = mne_read_events(eventname);
        catch
            error(me,mne_omit_first_line(lasterr));
        end
        fprintf(1,'Binary event file %s read\n',eventname);
    else
        %
        %   Text file
        %
        try
            events = load(eventname);
        catch
            error(me,mne_omit_first_line(lasterr));
        end
        if size(events,1) < 1
            error(me,'No data in the event file');
        end
        %
        %   Convert time to samples if sample number is negative
        %
        for p = 1:size(events,1)
            if events(p,1) < 0
                events(p,1) = events(p,2)*raw.info.sfreq;
            end
        end
        %
        %    Select the columns of interest (convert to integers)
        %
        events = int32(events(:,[1 3 4]));
        %
        %    New format?
        %
        if events(1,2) == 0 && events(1,3) == 0
            fprintf(1,'The text event file %s is in the new format\n',eventname);
            if events(1,1) ~= raw.first_samp
                error(me,'This new format event file is not compatible with the raw data');
            end
        else
            fprintf(1,'The text event file %s is in the old format\n',eventname);
            %
            %   Offset with first sample
            %
            events(:,1) = events(:,1) + raw.first_samp;
        end
    end
end
%
%    Select the desired events
%
count = 0;
selected = [];
for p = 1:size(events,1)
    if events(p,2) == 0 && events(p,3) == event
        count = count + 1;
        selected = [ selected p];
    end
end
if count > 0
    fprintf(1,'%d matching events found\n',count);
else
    error(me,'No desired events found.');
end

for p = 1:count
    %
    %       Read a data segment
    %
    event_samp = events(selected(p),1);
    from = event_samp + tmin*raw.info.sfreq;
    to   = event_samp + tmax*raw.info.sfreq;
    try
        if p == 1
            [ epoch ] = fiff_read_raw_segment(raw,from,to,picks);
            times = double([(int32(from)-int32(event_samp)):(int32(to)-int32(event_samp))])/raw.info.sfreq;
        else
            [ epoch ] = fiff_read_raw_segment(raw,from,to,picks);
        end
        data(p).epoch = epoch;
        data(p).event = event;
        data(p).tmin  = (double(from)-double(raw.first_samp))/raw.info.sfreq;
        data(p).tmax  = (double(to)-double(raw.first_samp))/raw.info.sfreq;
    catch
        fclose(raw.fid);
        error(me,'%s',mne_omit_first_line(lasterr));
    end
end
fprintf(1,'Read %d epochs, %d samples each.\n',count,length(data(1).epoch));

return;

end
