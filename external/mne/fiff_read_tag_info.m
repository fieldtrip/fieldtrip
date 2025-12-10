function [tag] = fiff_read_tag_info(fid)
%
% [fid,dir] = fiff_open(fname)
%
% Open a fif file and provide the directory of tags
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.2  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.1  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%

me='MNE:fiff_read_tag_info';

FIFFV_NEXT_SEQ=0;

data     = fread(fid,4,'int');
if ~isempty(data)
  tag.kind = data(1);
  tag.type = data(2);
  tag.size = data(3);
  tag.next = data(4);
else
  warning(me,'Invalid tag found');
  tag.kind = [];
  tag.type = [];
  tag.size = [];
  tag.next = [];
end

if tag.next == FIFFV_NEXT_SEQ
  fseek(fid,tag.size,'cof');
elseif tag.next > 0
  fseek(fid,tag.next,'bof');
end

