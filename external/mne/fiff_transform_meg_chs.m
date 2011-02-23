function [res, count] = fiff_transform_meg_chs(chs,trans)
%
% [res, count] = fiff_transform_meg_chs(chs,trans)
%
% Move to another coordinate system in MEG channel channel info
% Count gives the number of channels transformed
%
% NOTE: Only the coil_trans field is modified by this routine, not
% loc which remains to reflect the original data read from the fif file
%
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.7  2006/05/03 19:09:03  msh
%   Fixed even more compatibility issues.
%
%   Revision 1.6  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.5  2006/04/21 17:31:07  msh
%   Added the examples.
%   Modified the formatting of some informative output.
%
%   Revision 1.4  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.3  2006/04/17 15:01:34  msh
%   More small improvements.
%
%   Revision 1.2  2006/04/17 11:52:15  msh
%   Added coil definition stuff
%
%   Revision 1.1  2006/04/13 22:37:03  msh
%   Added head_head_trans field to info.
%
%   Revision 1.1  2006/04/13 21:20:06  msh
%   Added new routines for the channel information transformation.
%
%

me='MNE:fiff_transform_meg_chs';

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end

if nargin ~= 2
    error(me,'Wrong number of arguments');
 end

res = chs;
if isempty(trans)
   return;
end

count=0;
t   = trans.trans;
for k = 1:length(res)
    if res(k).kind == FIFF.FIFFV_MEG_CH || res(k).kind == FIFF.FIFFV_REF_MEG_CH
       if res(k).coord_frame == trans.from && ~isempty(res(k).coil_trans)
           res(k).coil_trans  = t*res(k).coil_trans;
           res(k).coord_frame = trans.to;
           count = count + 1;
       end
    end
end

if count > 0
    fprintf(1,'\t%d MEG channel locations transformed\n',count);
end

return;

end
