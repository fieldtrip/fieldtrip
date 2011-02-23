function fiff_write_raw_segment(fname, raw, from, to, sel, drop_small_buffer, buffer_size)
%   FIFF_WRITE_RAW_SEGMENT   Write chunck of raw data to disk
%       [] = FIFF_WRITE_RAW_SEGMENT(FNAME, RAW, FROM, TO, SEL)
%
%   The functions reads data from a file specified by raw
%   which is obtained with fiff_setup_read_raw
%
% fname                - the name of the file where to write
% raw                  - structure returned by fiff_setup_read_raw
% from                 - first sample to include. If omitted, defaults to the
%                        first sample in data
% to                   - last sample to include. If omitted, defaults to the last
%                        sample in data
% sel                  - optional channel selection vector
% drop_small_buffer    - optional bool to say if the last data buffer is dropped
%                        to make sure all buffers have the same size
%                        (required by maxfilter)
% buffer_size          - float (size of data buffers)

%
%   Author : Alexandre Gramfort, MGH Martinos Center
%   License : BSD 3 - clause
%

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end
%
me = 'MNE:fiff_write_raw_segment';
%
if nargin < 2
    error(me, 'Incorrect number of arguments');
end
if nargin < 3 | isempty(from)
    from = raw.first_samp;
end
if nargin < 4 | isempty(to)
    to = raw.last_samp;
end
if nargin < 5 | isempty(sel)
    sel = 1:raw.info.nchan;
end
if nargin < 6
    drop_small_buffer = false;
end
if nargin < 7
    buffer_size_sec = 25; % read by chunks of 30 seconds
    buffer_size = ceil(buffer_size_sec * raw.info.sfreq);
end
%
[outfid, cals] = fiff_start_writing_raw(fname, raw.info, sel);
%
first_buffer = true;
for first = from:buffer_size:to
    last = first + buffer_size - 1;
    if last > to
        last = to;
    end
    try
        [ data, times ] = fiff_read_raw_segment(raw, first, last, sel);
    catch
        fclose(raw.fid);
        fclose(outfid);
        error(me, '%s', mne_omit_first_line(lasterr));
    end
    if drop_small_buffer && first_buffer == false && length(times) < buffer_size
        fprintf(1, 'Skipping due to small buffer ... [done]\n');
        break
    end
    %
    %   You can add your own miracle here
    %
    fprintf(1, 'Writing...');
    if first_buffer
        if first > 0
            fiff_write_int(outfid, FIFF.FIFF_FIRST_SAMPLE, first);
        end
        first_buffer = false;
    end
    fiff_write_raw_buffer(outfid, data, cals);
    fprintf(1, '[done]\n');
end

fiff_finish_writing_raw(outfid);
