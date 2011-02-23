function [data,fid] = mne_read_epoch(epoch_info,which,prev_fid)
%
% [data,fid] = mne_read_epoch(epoch_info,which,prev_fid)
%
% Reads an epoch from a binary file produced by mne_epochs2mat
%
% epoch_info - The data structure read from the epoch data description file
% which      - Which epoch to read
% prev_fid   - Open file id from previous call
%              if prev_fid < 0 or missing, the file will be opened
%              The the current file id will be returned in the
%              output argument fid, if present. If this argument is
%              missing, file will be close upon exit from this function.
%
% The data will contain nchan x ntimes calibrated values
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%     $Id$
%     
%     Revision 1.6  2006/05/03 19:03:19  msh
%     Eliminated the use of cast function for Matlab 6.5 compatibility
%
%     Revision 1.5  2006/04/23 15:29:41  msh
%     Added MGH to the copyright
%
%     Revision 1.4  2006/04/10 23:26:54  msh
%     Added fiff reading routines
%
%     Revision 1.3  2006/02/21 03:35:35  msh
%     Improved help text of mne_read_epoch.m
%
%     Revision 1.2  2006/02/21 03:26:02  msh
%     Improved mne_read_epoch.m
%
%     Revision 1.1  2006/02/20 15:45:05  msh
%     Added mne_find_channel.m and mne_read_epoch.m
%
%
me='MNE:mne_read_epoch';
if(nargin ~= 2 && nargin ~= 3)
   error(me,'Usage : [data,fid] = mne_read_epoch(epoch_info,which,prev_fid)');
end

if (which < 0) 
   error(me,'Epoch number must be positive');
end

if (which > epoch_info.nepoch)
    error(me,'Epoch number too large');
 end

if (epoch_info.epochs(which,1) ~= 4)
   error(me,'Epoch is not in float format');
end

% open it as a big-endian file
if (nargin < 3 || prev_fid < 0) 
   [fid,message] = fopen(epoch_info.epoch_file, 'rb', 'b') ;
   if (fid < 0)
      error(me,message);
   end
else
   fid = prev_fid;
end

pos   = double(epoch_info.epochs(which,2));
nsamp = double(epoch_info.epochs(which,5));
nchan = double(epoch_info.nchan);

if (fseek(fid,pos,'bof') < 0)
  error(me,'Could not position file to the epoch');
end

data = fread(fid, [nchan,nsamp], 'float','ieee-be');
cal = diag(epoch_info.ch_cals(:,1))*diag(epoch_info.ch_cals(:,2));
data = cal*data;

if (nargout == 1)
   fclose(fid);
end

return;






