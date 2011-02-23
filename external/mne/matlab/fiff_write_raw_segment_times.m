function fiff_write_raw_segment_times(fname, raw, from, to, sel, drop_small_buffer, buffer_size_sec)
%   FIFF_WRITE_RAW_SEGMENT_TIMES   Write chunck of raw data to disk
%       [] = FIFF_WRITE_RAW_SEGMENT_TIMES(FNAME, RAW, FROM, TO, SEL)
%
%   The functions reads data from a file specified by raw
%   which is obtained with fiff_setup_read_raw
%
% fname              - the name of the file where to write
% raw                - structure returned by fiff_setup_read_raw
% from               - starting time of the segment in seconds
% to                 - end time of the segment in seconds
% sel                - optional channel selection vector
% drop_small_buffer  - optional bool to say if the last data buffer is dropped
%                      to make sure all buffers have the same size
%                      (required by maxfilter)
% buffer_size_sec    - float (size of data buffers in seconds)

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
    from = raw.first_samp / raw.info.sfreq;
end
if nargin < 4 | isempty(to)
    to = raw.last_samp / raw.info.sfreq;
end
if nargin < 5
    sel = 1:raw.info.nchan;
end
if nargin < 6
    drop_small_buffer = false;
end
if nargin < 7
    buffer_size_sec = 25;
end

buffer_size = ceil(buffer_size_sec * raw.info.sfreq);

%
%
%   Convert to samples
%
from = floor(from * raw.info.sfreq);
to   = ceil(to * raw.info.sfreq);

fiff_write_raw_segment(fname, raw, from, to, sel, drop_small_buffer, buffer_size);
