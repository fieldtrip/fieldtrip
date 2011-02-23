function mne_fwrite3(fid, val)
%
% mne_fwrite(fid, val)
% write a 3 byte integer to a file
% 

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

b1 = bitand(bitshift(val, -16), 255) ;
b2 = bitand(bitshift(val, -8), 255) ;
b3 = bitand(val, 255) ; 
fwrite(fid, b1, 'uchar') ;
fwrite(fid, b2, 'uchar') ;
fwrite(fid, b3, 'uchar') ;


