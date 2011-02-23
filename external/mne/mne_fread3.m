function [retval] = mne_fread3(fid)
%
% [retval] = mne_fread3(fid)
% read a 3 byte integer out of a file
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
% Revision 1.3  2006/04/23 15:29:40  msh
% Added MGH to the copyright
%
% Revision 1.2  2006/04/10 23:26:54  msh
% Added fiff reading routines
%
% Revision 1.1  2005/11/21 02:15:51  msh
% Added more routines
%
%

b1 = fread(fid, 1, 'uchar') ;
b2 = fread(fid, 1, 'uchar') ;
b3 = fread(fid, 1, 'uchar') ;
retval = bitshift(b1, 16) + bitshift(b2,8) + b3 ;

