function fiff_write_events(filename,eventlist,mappings)
%
% fiff_write_events(filename,eventlist,mappings)
%
% Write an event list into a fif file, and include an optional description
% of the event ids. This function has been adjusted by Jan Mathijs Schoffelen
% from mne_write_events, with the intention to make a writing function that
% is symmetric in its behavior w.r.t. fiff_read_events (which can read the
% mappings). The filename argument can be a string, or a file identifier to
% an open (for writing) fif-file.
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.2  2008/10/31 13:07:05  msh
%   Added mne_make_combined_event_file function
%
%   Revision 1.1  2008/06/16 17:27:50  msh
%   Added mne_read_events and mne_write_events functions
%

global FIFF;
if isempty(FIFF)
  FIFF = fiff_define_constants();
end

me = 'MNE:fiff_write_events';

if nargin < 3
  mappings = '';
end

eventlist = reshape(eventlist',numel(eventlist),1);
%
%   Start writing...
%
if ischar(filename)
  fid = fiff_start_file(filename);
else
  % assume the supplied filename to be a file identifier to an open file
  fid = filename;
end

fiff_start_block(fid,FIFF.FIFFB_MNE_EVENTS);
    
fiff_write_int(fid,FIFF.FIFF_MNE_EVENT_LIST,eventlist);
if ~isempty(mappings)
  % write a string that describes the mappings for numbers to a possibly
  % more informative string
  fiff_write_string(fid,FIFF.FIFF_DESCRIPTION,mappings);
end

fiff_end_block(fid,FIFF.FIFFB_MNE_EVENTS);

if ischar(filename)
  fiff_end_file(fid);
end
