function mne_make_combined_event_file(rawname,eventname,include,all,threshold)
%
%   mne_make_combined_event_file(rawname,eventname,include,all,threshold)
%
%   rawname     Name of the raw data file to scan
%   eventname   Name of the text format event file to output
%   include     Stimulus channel names to combine
%
%               This defaults to STI 001...STI 006
%
%   all         If true, include all trigger line transitions in the file
%               instead of the leading edges only
%   threshold   Threshold for detection of transition between inactive and active states
%
%   Create both a fif and eve format event file combining STI 001...STI 006
%   This function facilitates processing of Neuromag 122 data which do not
%   contain a composite trigger channel
%
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.4  2009/03/12 11:04:22  msh
%   Included a threshold parameter
%
%   Revision 1.3  2009/03/09 19:24:16  msh
%   Fixed errors in making the event files and allowed the text event file
%   not to be specified
%
%   Revision 1.2  2009/01/04 14:45:35  msh
%   The me string was wrong.
%
%   Revision 1.1  2008/10/31 13:07:04  msh
%   Added mne_make_combined_event_file function
%
%

%
me='MNE:mne_make_combined_event_file';
%
%   Setup for reading the raw data
%
try
    raw = fiff_setup_read_raw(rawname);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end
%
%   Set up pick list
%
if nargin < 3 || isempty(include)
    include{1} = 'STI 001';
    include{2} = 'STI 002';
    include{3} = 'STI 003';
    include{4} = 'STI 004';
    include{5} = 'STI 005';
    include{6} = 'STI 006';
end
if nargin < 4
    all = false;
end
if nargin < 5
    threshold = 0.1;
end
want_meg   = false;
want_eeg   = false;
want_stim  = false;
%
picks = fiff_pick_types(raw.info,want_meg,want_eeg,want_stim,include,raw.info.bads);

from = raw.first_samp;
to   = raw.last_samp;
%
%   Read a data segment
%   times output argument is optional
%
try
    [ data, times ] = fiff_read_raw_segment(raw,from,to,picks);
catch
    fclose(raw.fid);
    error(me,'%s',mne_omit_first_line(lasterr));
end
fprintf(1,'Read %d samples.\n',size(data,2));
%
%   Remember to close the file descriptor
%
fclose(raw.fid);
%
%   Make the combined channel
%
samples=[from:to];
comb=zeros(1,size(data,2));
for j = 1:size(data,1)
    for k = 1:size(data,2)
        if data(j,k) > threshold
            data(j,k) = 1;
            comb(k) = comb(k) + 2^(j-1);
        end
    end
end
%
%   Write the text file
%
if exist('eventname','var') && ~isempty(eventname)
    fd = fopen(eventname,'w');
    if fd < 0
        error(me,'Cannot open file %s',eventname);
    end
else
    fd = -1;
end

p = 0;
q = 0;
for k = 2:size(data,2)
    if comb(k) ~= comb(k-1)
        if fd >= 0 && (all || comb(k-1) == 0)
            fprintf(fd,'%6d %-10.4f %3d %3d\n',samples(k)-samples(1),times(k)-times(1),comb(k-1),comb(k));
            q = q + 1;
        end
        p = p + 1;
        events(p,1:3) = [ samples(k) comb(k-1) comb(k) ];
    end
end
if fd >= 0
    fprintf(1,'Wrote text event file %s (%d events)\n',eventname,q);
    fclose(fd);
end
%
%     Write the binary file
%
eventname_fif = strrep(rawname,'.fif','-eve.fif');
try
    mne_write_events(eventname_fif,events);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end
fprintf(1,'Wrote binary event file %s (%d events)\n',eventname_fif,size(events,1));

return;

end
