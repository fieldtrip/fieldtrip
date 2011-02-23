function fiff_write_name_list(fid,kind,data)
%
% fiff_write_name_list(fid,kind,mat)
% 
% Writes a colon-separated list of names
%
%     fid           An open fif file descriptor
%     kind          The tag kind
%     data          An array of names to create the list from
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.4  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.3  2006/04/12 10:29:02  msh
%   Made evoked data writing compatible with the structures returned in reading.
%
%   Revision 1.2  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%   Revision 1.1  2005/12/05 16:01:04  msh
%   Added an initial set of fiff writing routines.
%
%


me='MNE:fiff_write_name_list';

if nargin ~= 3
        error(me,'Incorrect number of arguments');
end

all=data{1};
for k = 2:size(data,2)
    all=strcat(all,':',data{k});
end
fiff_write_string(fid,kind,all);


