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

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end

FIFFV_NEXT_SEQ=0;

me='MNE:fiff_read_tag';

tag.kind = fread(fid,1,'int');
tag.type = fread(fid,1,'int');
tag.size = fread(fid,1,'int');
tag.next = fread(fid,1,'int');

if tag.next == FIFFV_NEXT_SEQ
  fseek(fid,tag.size,'cof');
elseif tag.next > 0
  fseek(fid,tag.next,'bof');
end

