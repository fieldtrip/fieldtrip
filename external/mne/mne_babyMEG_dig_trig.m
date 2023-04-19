function mne_baby_meg_dig_trig(infile,outfile,threshold,invert,include,want_eeg)
%
% function mne_baby_meg_dig_trig(infile,outfile,threshold,invert,want_eeg);
%
% Read and write raw data in 60-sec blocks
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%


global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end
%
me = 'MNE:mne_baby_meg_dig_trig';
%
if nargin < 2
    error(me,'Incorrect number of arguments');
end
if nargin < 3
   threshold = 4.8;
end
if nargin < 4
   invert = false;
end
if nargin < 5
   include = [];
end
if nargin < 6
    want_eeg = true;
end
%
%   Setup for reading the raw data
%
try
    raw = fiff_setup_read_raw(infile);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end
%
%raw.info.projs = [];
%
%   Set up pick list: MEG + STI 014 - bad channels
%
%
want_meg   = true;
want_stim  = false;
if isempty(include)
   include{1} = 'TRG001';
   include{2} = 'TRG002';
   include{3} = 'TRG003';
   include{4} = 'TRG004';
   include{5} = 'TRG005';
   include{6} = 'TRG006';
   include{7} = 'TRG007';
   include{8} = 'TRG008';
end
try
    picks = fiff_pick_types(raw.info,want_meg,want_eeg,want_stim,include,[]);
catch
    error(me,'%s (channel list may need modification)',mne_omit_first_line(lasterr));
end

fprintf(1,'Using threshold %.3f\n',threshold);
if invert
   fprintf(1,'Invert polarity\n');
else
   fprintf(1,'Keep polarity\n');
end
try
    dtrig = get_event_ch(infile,include,true,threshold,invert);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end
fprintf(1,'Digital trigger channel data ready.\n');
%
%
% Only one trigger channel remains and it will be called DTRG001
%
picksout = picks(1:end-length(include)+1);
raw.info.chs(picksout(end)).ch_name = 'DTRG001';
[outfid,cals] = fiff_start_writing_raw(outfile,raw.info,picksout);
%
%   Set up the reading parameters
%
from        = raw.first_samp;
to          = raw.last_samp;
quantum_sec = 10;
quantum     = ceil(quantum_sec*raw.info.sfreq);
%
%   To read the whole file at once set
%
%quantum     = to - from + 1;
%
%
%   Read and write all the data 
%   Replace the trigger channels with the single digital trigger channel
%
first_buffer = true;
for first = from:quantum:to
    last = first+quantum-1;
    if last > to
        last = to;
    end
    try
        [ data, times ] = fiff_read_raw_segment(raw,first,last,picks);
    catch
        fclose(raw.fid);
        fclose(outfid);
        error(me,'%s',mne_omit_first_line(lasterr));
    end
    %
    %   You can add your own miracle here
    %
    fprintf(1,'Writing...');
    if first_buffer
        first_first = first;
        if first > 0
            fiff_write_int(outfid,FIFF.FIFF_FIRST_SAMPLE,first);
        end
        first_buffer = false;
    end
    data(end-length(include)+1,:) = dtrig(first-first_first+1:last-first_first+1);
    fiff_write_raw_buffer(outfid,data(1:end-length(include)+1,:),cals);
    fprintf(1,'[done]\n');
end

fiff_finish_writing_raw(outfid);

return;


    function [comb] = get_event_ch(rawname,include,all,threshold,invert)
        %
        %   get_event_ch(rawname,include,all,threshold)
        %
        %   rawname     Name of the raw data file to scan
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
        me='MNE:get_event_ch';
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
        if nargin < 2 || isempty(include)
            include{1} = 'STI 001';
            include{2} = 'STI 002';
            include{3} = 'STI 003';
            include{4} = 'STI 004';
            include{5} = 'STI 005';
            include{6} = 'STI 006';
        end
        if nargin < 3
            all = false;
        end
        if nargin < 4
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
        fprintf(1,'Read %d samples of trigger inputs.\n',size(data,2));
        %
        %   Remember to close the file descriptor
        %
        %
        %   Make the combined channel
        %
        samples=[from:to];
        if invert
           comb=zeros(1,size(data,2));
           for j = 1:size(data,1)
              for k = 1:size(data,2)
                 if data(j,k) < threshold
                    data(j,k) = 1;
                    comb(k) = comb(k) + 2^(j-1);
                 end
              end
           end
        else
           comb=zeros(1,size(data,2));
           for j = 1:size(data,1)
              for k = 1:size(data,2)
                 if data(j,k) > threshold
                    data(j,k) = 1;
                    comb(k) = comb(k) + 2^(j-1);
                 end
              end
           end
        end
        
        return;

