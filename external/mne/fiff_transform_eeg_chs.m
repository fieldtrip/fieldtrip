function [res, count] = fiff_transform_eeg_chs(chs,trans)
%
% [res, count] = fiff_transform_eeg_chs(chs,trans)
%
% Move to another coordinate system in EEG channel channel info
% Count gives the number of channels transformed
%
% NOTE: Only the eeg_loc field is modified by this routine, not
% loc which remains to reflect the original data read from the fif file
%
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.7  2008/11/17 21:45:56  msh
%   Fixed error in transforming the EEG location data to another coordinate frame
%
%   Revision 1.6  2006/05/03 18:53:05  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.5  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.4  2006/04/21 17:31:07  msh
%   Added the examples.
%   Modified the formatting of some informative output.
%
%   Revision 1.3  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.2  2006/04/17 15:01:34  msh
%   More small improvements.
%
%   Revision 1.1  2006/04/17 11:52:15  msh
%   Added coil definition stuff
%
%

me='MNE:fiff_transform_eeg_chs';

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
%
%   Output unaugmented vectors from the transformation
%
t   = trans.trans(1:3,:);
for k = 1:length(res)
    if res(k).kind == FIFF.FIFFV_EEG_CH
       if res(k).coord_frame == trans.from && ~isempty(res(k).eeg_loc)
           %
           % Transform the augmented EEG location vectors
           %
           for p = 1:size(res(k).eeg_loc,2)
               res(k).eeg_loc(:,p) = t*[ res(k).eeg_loc(:,p) ; 1 ];
           end
           count = count + 1;
           res(k).coord_frame = trans.to;
       end
    end
end

if count > 0
    fprintf(1,'\t%d EEG electrode locations transformed\n',count);
end

return;

end
